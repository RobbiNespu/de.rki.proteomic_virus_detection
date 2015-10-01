# General Information

This is the Repo of the source for nodes of a KNIME workflow to detect virusses in human genome data.

It was part of a student project I did December 2014 to March 2015 in the Bioinformatics group of Bernhard Renard at the Robert-Koch-Institute in Berlin, Germany.

The aim of the project was to detect virus families and later strains in MS2 data of infected human cells.

# Build Process

The source is build using [GenericKnimeNodes](https://github.com/genericworkflownodes/GenericKnimeNodes) into a KNIME plugin that can then be used to build workflows. A building instruction can be found in the [Seqan Documentation](http://seqan.readthedocs.org/en/latest/HowTo/GenerateKnimeNodesExternalTools.html).

Note that the folders in the plugin_source directory have to be zipped to be used (eg. payload/binaries_win_32.zip). For Windows there are two versions of the plugin possible: a normal python version and a version of precompiled scritps that can be created with the setup_py2exe.py script. Second version is much larger in size but does not have python as a dependency which has to be set specifically by the enduser under payload/binaries.ini.