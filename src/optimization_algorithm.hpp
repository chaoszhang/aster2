#ifndef OPTIMIZATION_ALGORITHM_HPP
#define OPTIMIZATION_ALGORITHM_HPP

#include "common.hpp"
#include "placement_algorithm.hpp"

namespace optimization_algorithm {

	using std::size_t;
	using std::vector;

	template<placement_algorithm::PLACEMENT_ALGORITHM PlacementAlgorithm> class RecursivePlacement : public common::LogInfo{
	public:
		RecursivePlacement(int verbose = common::LogInfo::DEFAULT_VERBOSE) : LogInfo(verbose){}

		common::AnnotatedBinaryTree& operator()(typename PlacementAlgorithm::Color& color, typename PlacementAlgorithm::ThreadPool& threadPool, common::AnnotatedBinaryTree& tree, std::ranges::range auto const& taxa) {
			size_t nLeaves = tree.leaves().size();
			LogInfo const& vlog = log.vlog<1>();
			for (size_t iTaxon : taxa) {
				if (nLeaves < 3) tree.emplaceRoot(common::AnnotatedBinaryTree::LEAF_ID, iTaxon);
				else {
					log << "Placing " << nLeaves + 1 << "-th taxon (" << common::taxonName2ID[iTaxon] << ")..." << std::endl;
					tree.displaySimpleNewick<double, double>(vlog << "Current tree: ", "", "") << std::endl;
					PlacementAlgorithm placement(color, threadPool, tree, verbose + 1);
					placement.place(iTaxon);
				}
				nLeaves++;
			}
			return tree;
		}
	};

	template<stepwise_colorable::STEPWISE_COLORABLE SC> class SimpleStepwiseColor : public common::LogInfo {
		using StepwiseColor = SC;
	public:
		SimpleStepwiseColor(int verbose = common::LogInfo::DEFAULT_VERBOSE) : LogInfo(verbose){}

		template<class Generator> common::AnnotatedBinaryTree operator()(typename StepwiseColor::SharedConstData const& data, size_t nTaxa, common::Random<Generator> &random, size_t nThreads = 1) {
			vector<size_t> taxa = random.randomTaxonOrder(nTaxa);
			StepwiseColor stepwiseColor(&data);
			thread_pool::ThreadPool<typename StepwiseColor::score_t, thread_pool::SimpleScheduler> tp(nThreads, 0, data.nElements);
			common::AnnotatedBinaryTree tree;
			optimization_algorithm::RecursivePlacement<placement_algorithm::StepwiseColorPlacement<placement_algorithm::StepwiseColorPlacementDefaultAttributes<StepwiseColor> > > rp(verbose);
			rp(stepwiseColor, tp, tree, taxa);
			tree.displaySimpleNewick<double, double>(log << "Final tree: ", "", "") << std::endl;
			placement_algorithm::StepwiseColorPlacement<placement_algorithm::StepwiseColorPlacementDefaultAttributes<StepwiseColor> > placement(stepwiseColor, tp, tree, verbose);
			log << "Final score: " << placement.scoreTree() << std::endl;
			return tree;
		}
	};
};

#endif