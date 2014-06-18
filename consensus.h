#include <seqan/sequence.h>
#include <seqan/store.h>
#include <seqan/graph_align.h>
#include <seqan/align.h>
#include <seqan/graph_types.h>
#include <seqan/graph_msa.h>
#include <seqan/realign.h>

#include "ag_to_frag_store.h"

using namespace seqan;

typedef String<Fragment<> > TFragments;

// Compute all-to-all unbanded overlap alignments of seqs.
template<typename TFragments, typename TSeq>
void computeOverlapAlignments(TFragments & frags,
                                Graph<Undirected<double> > & distances,
                                String<int> & scores,
                                StringSet<TSeq> & seqs)
{
    Score<int, Simple> scoringScheme(1, -1, -1);
    
    for (unsigned i = 0; i < length(seqs); ++i)
    {
        for (unsigned j = i + 1; j < length(seqs); ++j)
        {
            AlignConfig<true, true, true, true> alignConfig;

            StringSet<Dna5String, Dependent<> > pairSet;
            assignValueById(pairSet, seqs, i);
            assignValueById(pairSet, seqs, j);

            TFragments matches;
            int score = globalAlignment(matches, pairSet, scoringScheme, alignConfig);
            if (empty(matches)) // no overlap
                continue;

            // Compute alignment statistics.
            unsigned matchLen = 0;
            unsigned overlapLen = 0;
            unsigned alignLen = 0;
            getAlignmentStatistics(matches, pairSet, 0u, (unsigned)length(matches), matchLen, overlapLen, alignLen);

            if (matchLen < 13)
                continue;

            int quality = 100.0 * matchLen / overlapLen; // percentage of matches in overlap
            if (quality < 95)
                continue;

            // Append alignment and score.
            append(frags, matches);
            resize(scores, length(frags), score); // each fragment is scored with score

            // Add edge to distance matrix.
            addEdge(distances, i, j, quality);
        }
    }
}

template<typename TFragments, typename TSeq>
int multiReadAlignment(FragmentStore<void, MyStoreConfig_> & store,
                       TFragments & frags,
                       Graph<Undirected<double> > & distances,
                       String<int> & scores,
                       StringSet<TSeq> & seqs)
{
    typedef StringSet<Dna5String, Dependent<> > TDepReadSet;
    typedef Graph<Alignment<TDepReadSet, unsigned> > TInGraph;
    
    // Build the alignment graph.
    Score<int, Simple> msaScoringScheme(2, -6, -4, -9);
    TDepReadSet depSeqs(seqs);
    TInGraph inGraph(depSeqs);
    buildAlignmentGraph(frags, scores, inGraph, msaScoringScheme, ReScore());

    // Perform triplet library extension.
    tripletLibraryExtension(inGraph);

    // Compute guide tree.
    Graph<Tree<double> > guideTree;
    Graph<Undirected<double> > dCopy(distances);
    upgmaTree(dCopy, guideTree);

    // Perform progressive alignment.
    TInGraph graph(depSeqs);
    assignStringSet(graph, stringSet(inGraph));
    progressiveAlignment(inGraph, guideTree, graph);

    // Build resulting alignment graph.
    String<unsigned> component;
    String<unsigned> order;
    std::map<unsigned, unsigned> componentLength;
    SEQAN_ASSERT_NOT(empty(graph));
    if (!convertAlignment(graph, component, order, componentLength))
    {
        std::cerr << "ERROR: Built alignmetn graph was invalid.\n";
        return 1;
    }

    // Convert graph to FragmentStore.
    unsigned numComponents = length(order);
    alignmentGraphToFragmentStore(store, graph, distances, component, order, numComponents);
    
    return 0;
}

template <typename TFragmentStore>
void printStore(std::ostream & out, TFragmentStore const & storeC, unsigned contigID)
{
    TFragmentStore store(storeC);

    AlignedReadLayout layout;
    layoutAlignment(layout, store);
    for (unsigned i = 0; i < length(layout.contigRows); ++i)
    {
        if (contigID != (unsigned)-1 && contigID != i)
            continue;

        // Get coordinates to plot.
        __int64 l = 0;
        __int64 r = l;
        for (unsigned j = 0; j < length(layout.contigRows[i]); ++j)
        {
            unsigned id = back(layout.contigRows[i][j]);
            if (r < store.alignedReadStore[id].beginPos)
                r = store.alignedReadStore[id].beginPos;
            if (r < store.alignedReadStore[id].endPos)
                r = store.alignedReadStore[id].endPos;
        }

        out << ">multi-read-alignment: contig_" << i << "\n";
        printAlignment(out, Raw(), layout, store, i, l, r, 0, 1000);
    }
}

template<typename TSeq1, typename TSeq2>
int compute_consensus(StringSet<TSeq1> & consensusSeqs, StringSet<TSeq2> & seqs)
{
    // Compute overlap alignments from seqs, all-to-all and unbanded here.
    TFragments frags;
    String<int> scores;
    Graph<Undirected<double> > distances;
    _resizeWithRespectToDistance(distances, length(seqs));
    computeOverlapAlignments(frags, distances, scores, seqs);
    
    // Compute a multiple read alignment and write it to a FragmentStore.
    FragmentStore<void, MyStoreConfig_> store;
    if (multiReadAlignment(store, frags, distances, scores, seqs) == 1)
        return 1;

    // Realign FragmentStore, will also generate consensus sequences.
    for (unsigned i = 0; i < length(store.contigStore); ++i)
        reAlignment(store, i, 1, 10, false);

// **************************************************************
// ***********************  DEBUG CODE **************************

//    // Print realigned MSA with consensus.
//    std::cerr << "Re-aligned MSAs\n";
//    for (unsigned i = 0; i < length(store.contigStore); ++i)
//        printStore(std::cerr, store, i);

// **************************************************************

    // Collect the consensus sequences.
    for (unsigned i = 0; i < length(store.contigStore); ++i)
        appendValue(consensusSeqs, store.contigStore[i].seq);

    return 0;
}
