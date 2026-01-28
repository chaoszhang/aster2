#include "optimization_algorithm.hpp"

#ifdef CASTER
#include "caster.hpp"
namespace my_tool = caster;
#else
#include "caster.hpp"
namespace my_tool = caster;
#endif

using std::cerr;
using std::endl;
using namespace std::string_literals;

int main(int argc, char* argv[]) {
	using std::string;
	using std::size_t;

	std::pair<string, string> programNames = my_tool::Driver::programNames();
	ARG.set("SHORT_NAME", programNames.first); 
	ARG.set("FULL_NAME", programNames.second);

	ARG.addArgument('h', "help", "flag", "Display help message", 6, true);
	ARG.addArgument('i', "input", "string", "Input file path", 5);
	ARG.addArgument('o', "output", "string", "Output file path, print to stdout if not provided", 5, true);
	ARG.addArgument('t', "thread", "integer", "Number of threads", 4, true, true, "1");
	ARG.addArgument('a', "mapping", "string", "Mapping file path, a list of gene/speicesman name to taxon name maps, each line contains one gene/speicesman name followed by one taxon name separated by a space or tab", 3, true, false);
	ARG.addArgument('\0', "verbose", "integer", "Verbose level", 0, true, true, "2");
	ARG.addArgument('\0', "no-log", "flag", "Don't generate log file", 1, true);
	ARG.addArgument('\0', "log", "string", "Log file path", 0, true, true, "log.txt");
	ARG.addArgument('\0', "log-verbose", "integer", "Verbose level in log file", 0, true, true, "3");
	my_tool::Driver::addArguments();

	ARG.parse(argc, argv);
	common::LogInfo::setVerbose(ARG.has("no-log") ? nullptr : new std::ofstream(ARG.get<string>("log")), ARG.get<size_t>("log-verbose"), ARG.get<size_t>("verbose"));
	ARG.print();

	size_t nThreads = ARG.get<size_t>("thread");
	ARG << "#threads: " << nThreads << endl;

	auto stepwiseColorSharedConstData = std::visit([](auto&& arg) -> decltype(auto) {
		return std::forward<decltype(arg)>(arg);
	}, my_tool::Driver::getStepwiseColorSharedConstData());
	size_t nTaxa = common::taxonName2ID.nTaxa();
	ARG << "#taxa: " << nTaxa << endl;

	optimization_algorithm::Procedure<optimization_algorithm::DefaultProcedureAttributes<typename decltype(stepwiseColorSharedConstData)::ParentClass> > alg;
	std::shared_ptr<common::AnnotatedBinaryTree> tree = alg.constrainedDPProcedure(stepwiseColorSharedConstData, nTaxa, 4, nThreads, 0);
	
	tree->displaySimpleNewick<double, double>(ARG << "Final tree: ", "", "") << endl;


	if (ARG.has("output")) {
		std::ofstream fout(ARG.get<string>("output"));
		tree->displaySimpleNewick<double, double>(fout, "support", "length");
		fout << endl;
	}
	else {
		tree->displaySimpleNewick<double, double>(std::cout, "support", "length");
		std::cout << endl;
	}
	return 0;
}