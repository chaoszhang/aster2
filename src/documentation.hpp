#ifndef DOCUMENTATION_HPP
#define DOCUMENTATION_HPP

class DocumentationBase {
protected:
	virtual string introduction() const noexcept = 0;

	virtual string documentations() const noexcept {
		return R"YOHANETYO(# Documentations
- The rest of this documentation file
- Forums (feel free to ask questions or ask for help running ASTER2):
  - [User group discussions](https://groups.google.com/forum/#!forum/aster-users)
  - [ASTER2 issues page](https://github.com/chaoszhang/aster2/issues)
  - QQ group: 130635706

## Bug Reports
Contact ``chaozhang@pku.edu.cn``, [``aster-users@googlegroups.com``](https://groups.google.com/forum/#!forum/aster-users), or post on [ASTER2 issues page](https://github.com/chaoszhang/aster2/issues).
)YOHANETYO";
	}

	virtual string installation() const noexcept {
		return R"YOHANETYO(# INSTALLATION
```
make
```
)YOHANETYO";
	}

	virtual string input() const noexcept = 0;

	virtual string mapping() const noexcept {
		return R"YOHANETYO(## Mapping file
ASTER2 assumes that taxon names in the input are consistent throughout.
Otherwise, it will treat them as different taxa, unless you provide a mapping file.
The mapping file should be a text file where each line contains an alias followed by a taxon name, separated by a space or tab. For example:
```
alias_A1 taxon_A
alias_A2 taxon_A
alias_B1 taxon_B
alias_B2 taxon_B
alias_B3 taxon_B
...
```
ASTER2 will internally replace all occurrences of aliases with corresponding taxon names.
If an alias in the mapping file does not appear in the input, it will be ignored.
A mapping file is also useful when you don't want to use the original taxon names in the input as taxon names in the output.
A mapping file is convient when you have multiple speciesmen per taxon and you want to specify which speciesmen belong to which taxa.
)YOHANETYO";
	}

	virtual string output() const noexcept {
		return R"YOHANETYO(# OUTPUT
The output in is Newick format and gives:
* the species tree topology
* branch supports measured as local block bootstrap supports by default (>95.0 means good)
* It can also annotate branches with other quantities, such as quartet scores and bootstraps for all three topologies.
)YOHANETYO";
	}

	virtual string programName() const noexcept = 0;

	virtual string exampleInput() const noexcept = 0;

	virtual string execution() const noexcept {
		return std::format(R"YOHANETYO(# EXECUTION
ASTER2 currently has no GUI. You need to run it through the command-line. In a terminal/PowerShell, go to the directory (location) where you have downloaded ASTER2 and issue the following command:

```
bin/{0}
```

This will give you a list of options available. If you are using Windows, please replace `bin/{0}` with `.\exe\{0}.exe`.

To find the species tree with input from in a file called `INPUT_FILE`, use:
```
bin/{0} -i INPUT_FILE
```

For example if you want to run `{0}` with input `example/{1}`, then run

```
bin/{0} -i example/{1}
```

The results will be outputted to the standard output. To save the results in a file use the `-o OUTPUT_FILE` option:
```
bin/{0} -i INPUT_FILE -o OUTPUT_FILE
```

```
bin/{0} -i example/{1} -o example/{1}.stree
```

When using a mapping file, add `-a MAPPING_FILE`:
```
bin/{0} -i INPUT_FILE -o OUTPUT_FILE 
```

ASTER2 supports multi-threading. To run program with 4 threads, add `-t 4`:

```
bin/{0} -t 4 -i INPUT_FILE -o OUTPUT_FILE
```

ASTER2 has very good parrallel efficiency up to 64 cores when input data is large. In fact, it often experiences super-linear speedup with 16 cores or more. So feel free to use as many cores as you want.

ASTER2 also allows rooting at an given outgroup taxon:

```
bin/{0} --root YOUR_OUTGROUP INPUT_FILE
```
)YOHANETYO", programName(), exampleInput());
	}

	virtual string changes() const noexcept {
		std::ostringstream out;
		out << "# CHANGE LOG\n";
		ChangeLog::displayAll<true>(out, ChangeLog::shared.all);
		return out.str();
	}

public:
	virtual string operator()() const noexcept {
		return introduction() + "\n" + documentations() + "\n" + installation() + "\n" + input() + "\n" + mapping() + "\n" + output() + "\n" + execution() + "\n" + changes();
	}
};

#endif