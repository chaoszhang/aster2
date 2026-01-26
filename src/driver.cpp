#include "caster.hpp"
#include "placement_algorithm.hpp"
#include "optimization_algorithm.hpp"

using namespace caster;
using namespace placement_algorithm;
using namespace optimization_algorithm;

using std::cerr;
using std::endl;


int main(int argc, char* argv[]) {

	ChangeLog::displayAll<true>(std::cout, ChangeLog::shared.all);
	ChangeLog::displayAll<false>(std::cout, ChangeLog::shared.byCode["ChangeLog"]);
	using PlacementAttributes = StepwiseColorPlacementDefaultAttributes<StepwiseColor<StepwiseColorDefaultAttributes> >;

	StepwiseColor<StepwiseColorDefaultAttributes>::SharedConstData stepwiseColorSharedConstData;
	common::Random<std::mt19937> random;
	SimpleStepwiseColor<StepwiseColor<StepwiseColorDefaultAttributes> > simpleStepwiseColor(1);
	common::AnnotatedBinaryTree tree = simpleStepwiseColor(stepwiseColorSharedConstData, 3, random, 2);
	tree.verbose = 1;
	common::LogInfo::setVerbose(nullptr, 2, 2);

	tree.logging("verbose 0");
	tree.logging(endl);
	tree << "verbose 0" << endl;
	tree.logging<1>("verbose 1");
	tree.logging<1>(endl);
	return 0;
}