#ifndef AG_TO_FRAG_STORE_H_


struct MyStoreConfig_
{
typedef seqan::String<seqan::Dna5Q>	TReadSeq;
typedef seqan::String<seqan::Dna5Q>	TContigSeq;

typedef double	TMean;
typedef double	TStd;
typedef signed char	TMappingQuality;

typedef void	TReadStoreElementSpec;
typedef seqan::Owner<> TReadSeqStoreSpec;
typedef void	TMatePairStoreElementSpec;
typedef void	TLibraryStoreElementSpec;
typedef void	TContigStoreElementSpec;
typedef void	TContigFileSpec;
typedef void	TAlignedReadStoreElementSpec;
typedef seqan::Owner<>	TAlignedReadTagStoreSpec;
typedef void	TAnnotationStoreElementSpec;
    
    typedef seqan::Alloc<>	TReadNameSpec;
typedef seqan::Owner<>	TReadNameStoreSpec;
};

namespace seqan {

// ---------------------------------------------------------------------------
// Function alignmentGraphToSmoothFragmentStore()
// ---------------------------------------------------------------------------

template <typename TFragmentStore, typename TSequence, typename TCargo, typename TSetSpec, typename TSpec>
bool alignmentGraphToFragmentStore(TFragmentStore & store,
                                   seqan::Graph<seqan::Alignment<seqan::StringSet<TSequence, TSetSpec>, TCargo, TSpec> > const & g,
                                   seqan::Graph<seqan::Undirected<double> > const & distances,
                                   seqan::String<unsigned> const & component,
                                   seqan::String<unsigned> const & order,
                                   unsigned numComponents)
{
    // NOTE: seqToCluster is indexed by POSITION in the read set of g and not by the ID.

    using namespace seqan;

    typedef Graph<Alignment<StringSet<TSequence, TSetSpec>, TCargo, TSpec> > TAlignmentGraph;

    // Allocate information for which sequence supports the profile at which position.
    resize(store.readSeqStore, length(stringSet(g)));
    resize(store.readStore, length(stringSet(g)));
    resize(store.alignedReadStore, length(stringSet(g)));

    // We can fill the mate pair store here since we know that we have an even number of reads in the read set and the
    // read with a given id is part of the pair (id / 2) and is the (id % 2)-th read in the pair.
    resize(store.matePairStore, length(store.readStore));
    for (unsigned i = 0; i < length(store.matePairStore); ++i)
        store.matePairStore[i].readId[0] = i;

    // -----------------------------------------------------------------------
    // Get connected components of distances / read alignment clusters.
    // -----------------------------------------------------------------------

    // Each cluster corresponds to a contig.

    // A cluster is a CC in the graph where each sequences is a vertex and two vertices are connected if they have an
    // overlap alignment.
    String<unsigned> seqToCluster;
    unsigned numClusters = connectedComponents(distances, seqToCluster);
    resize(store.contigStore, numClusters);
    String<unsigned> contigLengths;
    resize(contigLengths, numClusters, 0);

    for (unsigned i = 0; i < numClusters; ++i)
    {
        std::stringstream ss;
        ss << "contig_" << i;
        appendValue(store.contigNameStore, ss.str());
    }

    // -----------------------------------------------------------------------
    // Visit components in topological order and generate profile sequences.
    // -----------------------------------------------------------------------

    // Get mapping from component to vertices.
    String<String<unsigned> > componentVertices;
    resize(componentVertices, numComponents);
    typedef typename Iterator<TAlignmentGraph, VertexIterator>::Type TVertexIterator;
    for (TVertexIterator itV(g); !atEnd(itV); goNext(itV))
        appendValue(componentVertices[getProperty(component, *itV)], *itV);

    // For each cluster, the currently overlapping reads.
    std::vector<std::set<unsigned> > activeReads(numClusters);

    std::vector<unsigned> gapCount(length(stringSet(g)), 0);

    // Iterate vertices in topological order.
    for (typename Iterator<String<unsigned> const, Rooted>::Type it = begin(order, Rooted()); !atEnd(it); goNext(it))
    {
        unsigned c = *it; // Current component.
        unsigned fLen = fragmentLength(g, front(componentVertices[c]));
        for (unsigned i = 1; i < length(componentVertices[c]); ++i)
            SEQAN_ASSERT_EQ(fragmentLength(g, front(componentVertices[c][0])),
                            fragmentLength(g, front(componentVertices[c][i])));
        unsigned cl = seqToCluster[idToPosition(stringSet(g), sequenceId(g, front(componentVertices[c])))]; // Current cluster/contig.

        // Update contig lengths.
        unsigned from = contigLengths[cl];
        contigLengths[cl] += fLen;

        // The currently active reads that we see in this round. Required for inserting gaps below.
        std::set<unsigned> seen;
        std::set<unsigned> done;

        // Insert gaps.
        typedef typename Iterator<String<unsigned>, Rooted>::Type TDescIt;
        for (TDescIt itV = begin(componentVertices[c], Rooted()); !atEnd(itV); goNext(itV))
        {
            unsigned idx = idToPosition(stringSet(g), sequenceId(g, *itV));
            seen.insert(idx);
            unsigned fBeg = fragmentBegin(g, *itV);

            // Register sequence as supporting in profile cl starting at position from in profile.
            if (fBeg == 0u)
            {
                store.readSeqStore[idx] = getValueById(stringSet(g), sequenceId(g, *itV));
                store.readStore[idx].matePairId = idx;
                SEQAN_ASSERT_NOT(empty(store.readSeqStore[idx]));
                activeReads[cl].insert(idx);
                store.alignedReadStore[idx].id = idx;
                store.alignedReadStore[idx].readId = idx;
                store.alignedReadStore[idx].contigId = cl;
                store.alignedReadStore[idx].beginPos = from;
                store.alignedReadStore[idx].endPos = from;
                store.alignedReadStore[idx].pairMatchId = idx / 2;
            }
            store.alignedReadStore[idx].endPos = from + fLen;

            unsigned fEnd = fBeg + fLen;

            if (fEnd == length(stringSet(g)[idx]))
                done.insert(idx);
        }

        // Get not seen reads.
        typedef std::set<unsigned>::iterator TSetIt;
        std::set<unsigned> notSeen;
        for (TSetIt it = activeReads[cl].begin(); it != activeReads[cl].end(); ++it)
            notSeen.insert(*it);
        for (TSetIt it = seen.begin(); it != seen.end(); ++it)
            notSeen.erase(*it);
        // Insert gaps into these reads.
        for (TSetIt itS = notSeen.begin(); itS != notSeen.end(); ++itS)
        {
            typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
            typedef typename Value<TAlignedReadStore>::Type TAlignedRead;
            typedef typename TAlignedRead::TGapAnchors TGapAnchors;
            typedef typename TFragmentStore::TReadSeq TReadSeq;
            SEQAN_ASSERT_NOT(empty(store.readSeqStore[*itS]));
            Gaps<TReadSeq, AnchorGaps<TGapAnchors> > gaps(store.readSeqStore[*itS], store.alignedReadStore[*itS].gaps);
            insertGaps(gaps, from - store.alignedReadStore[*itS].beginPos, fLen);
            store.alignedReadStore[*itS].endPos += fLen;
            gapCount[*itS] += fLen;
        }

        // Deactive done reads.
        for (TSetIt it = done.begin(); it != done.end(); ++it)
            activeReads[cl].erase(*it);
    }
 
// #if SEQAN_ENABLE_DEBUG
    {
        // Check for consistency.
        typedef typename TFragmentStore::TAlignedReadStore TAlignedReadStore;
        typedef typename Iterator<TAlignedReadStore, Standard>::Type TAlignedReadIter;
        typedef typename TFragmentStore::TReadSeq TReadSeq;

        TAlignedReadIter itEnd = end(store.alignedReadStore, Standard());
        for (TAlignedReadIter it2 = begin(store.alignedReadStore, Standard()); it2 != itEnd; ++it2)
        {
            typedef Gaps<TReadSeq, AnchorGaps<String<typename TFragmentStore::TReadGapAnchor> > > TReadGaps;
            TReadGaps readGaps(store.readSeqStore[it2->readId], it2->gaps);
            SEQAN_ASSERT_EQ(length(readGaps) - length(store.readSeqStore[it2->readId]), gapCount[it2->readId]);
            if ((unsigned)abs(it2->endPos - it2->beginPos) != length(readGaps))
            {
                SEQAN_FAIL("Inconsistent begin/endPos");
            }
        }
    }
// #endif // #if SEQAN_ENABLE_DEBUG

    return true;
}

} // namespace seqan

#endif // #ifndef AG_TO_FRAG_STORE_H_
