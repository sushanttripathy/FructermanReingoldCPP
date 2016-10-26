#include <iostream>
#include "fructermanreingold.h"
#include "csv.h"
#include "liboptions.hpp"

int OptionsHelp() {
    std::cout
            << "Usage : FR <options>                                            " << std::endl
            << std::endl
            << "Options :                                                       " << std::endl
            << std::endl
            << "--edge-list         <edge list csv path>                        " << std::endl
            << "--output-coords     <nodes coordinates output csv path>         " << std::endl
            << "--node-weights      <nodes weights csv path (optional)>         " << std::endl
            << std::endl
            << "--width             <int default 100>                           " << std::endl
            << "--height            <int default 100>                           " << std::endl
            << "--iterations        <int default 300>                           " << std::endl
            << "--initial-temp      <float (optional)>                          " << std::endl
            << "--temp-decay-const  <float default 0.9>                         " << std::endl
            << std::endl
            << "--ignore-edge-weights                                           " << std::endl
            << std::endl
            << "--verbosity-level   <int>                                       " << std::endl
            << std::endl;

    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        OptionsHelp();
    } else {
        Options o(argc, argv, "--");
        if (o.isSet("--edge-list") && o.isSet("--output-coords")) {

            int width = o.GetInteger("--width", 100);
            int height = o.GetInteger("--height", 100);


            FructermanReingold F(width, height);

            try {
                io::CSVReader<3> in_csv(o.GetString("--edge-list"));
                in_csv.read_header(io::ignore_extra_column | io::ignore_missing_column | io::ignore_no_column, "source",
                                   "target", "weight");

                if (!in_csv.has_column("source") || !in_csv.has_column("target")) {
                    std::cerr << " This software expects the edge list csv to contain these columns at minimum: "
                              << std::endl;
                    std::cerr << " source, target" << std::endl;
                    std::cerr << " It can optionally contain a 'weight' column." << std::endl;
                    exit(0);
                }

                int64_t source, target;
                float weight = 1.0;
                while (in_csv.read_row(source, target, weight)) {
                    F.addEdge(source, target, weight);
                }

            } catch (std::exception &e) {
                std::cerr << e.what() << std::endl;

                exit(0);
            }

            if (o.isSet("--node-weights")) {
                try {
                    io::CSVReader<2> in_csv(o.GetString("--node-weights"));
                    in_csv.read_header(io::ignore_extra_column | io::ignore_missing_column | io::ignore_no_column,
                                       "node id", "weight");

                    if (!in_csv.has_column("node id") || !in_csv.has_column("weight")) {
                        std::cerr << " This software expects the node weights csv to contain these columns at minimum: "
                                  << std::endl;
                        std::cerr << " node id, weight" << std::endl;
                    }

                    int64_t node_id;
                    float weight = 1.0;
                    while (in_csv.read_row(node_id, weight)) {
                        F.addNodeWeight(node_id, weight);
                    }

                } catch (std::exception &e) {
                    std::cerr << e.what() << std::endl;
                    exit(0);
                }
            }

            try {
                long num_iterations = o.GetInteger("--iterations", 300);
                float init_temp = o.GetReal("--initial-temp", -1.0);
                float decay_constant = o.GetReal("--temp-decay-const", 0.9);

                if (o.isSet("--ignore-edge-weights")) {
                    F.setIgnoreEdgeWeights(true);
                }

                if (o.isSet("--verbosity-level")) {
                    F.setVerbose(o.GetInteger("--verbosity-level", 0));
                }

                F.run(num_iterations, init_temp, decay_constant);
                F.saveCSV(o.GetString("--output-coords"));
            } catch (std::exception &e) {
                std::cerr << e.what() << std::endl;
                exit(0);
            }


        } else if (!(o.isSet("--edge-list") || o.isSet("--output-coords"))) {
            std::cerr << " Both the edge list csv and output coordinates need to be specified." << std::endl;
            std::cerr
                    << " Please sepcify them as such: --edge-list /path/to/file1name.csv --output-coords /path/to/file2name.csv"
                    << std::endl;
        } else if (!o.isSet("--edge-list")) {
            std::cerr << " Please specify the edge list csv using  : --edge-list /path/to/filename.csv" << std::endl;
        } else {
            std::cerr << " Please specify the output coordinates csv using  : --output-coords /path/to/filename.csv"
                      << std::endl;
        }

    }
    return 0;
}