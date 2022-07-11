#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>

#include <boost/range/irange.hpp>
#include <boost/program_options.hpp> 

#include "common.h"
#include "total_complex.h"
#include "utils.h"
#include "dpc.h"

bool parseOptions(int argc, char** argv, std::string& posfname, 
    Integer& max_dim, Integer& min_time, Integer& max_time, DynamicPC& dpc) {

    try { 
        std::string intervals_str;

        namespace po = boost::program_options; 
        po::options_description shown_desc("Options"); 
        shown_desc.add_options() 
            (",d", po::value<Integer>(&max_dim), 
            "Maximum dimension for Rips simplices\n")
            (",s", po::value<Integer>(&min_time), 
            "Starting time considered for the DPC\n")
            (",e", po::value<Integer>(&max_time), 
            "Ending time considered for the DPC\n")
            ("fzz", "Compute barcodes using FastZigzag for timing comparison\n")
            ("eop,E", "Print edge operations\n")
            ("sop", "Print simplex operations (automatically turns on '-E')\n")
            ("help,h", "Print help messages\n");

        po::options_description hidden_desc("Hidden options"); 
        hidden_desc.add_options() 
            ("filename", po::value<std::string>(&posfname)->required(), "filename");

        po::options_description desc("All options"); 
        desc.add(shown_desc).add(hidden_desc);

        po::positional_options_description positional_options; 
        positional_options.add("filename", 1);

        po::variables_map vm;

        try { 
            po::store(po::command_line_parser(argc, argv).
                options(desc).positional(positional_options).run(), vm);

            if (vm.count("help")) { 
                cout << endl << "USAGE: ./dpc_vine [OPTIONS] [FILENAME]" << endl << endl
                    << shown_desc; 
                exit(0); 
            }

            po::notify(vm);

        } catch(boost::program_options::required_option& e) {  

            cout << "ERROR: " << e.what() << endl << endl; 
            return false; 

        } catch(boost::program_options::error& e) { 

            cout << "ERROR: " << e.what() << endl << endl; 
            return false; 

        }

        if (vm.count("eop")) { dpc.print_e_filt_op = true; }
        if (vm.count("sop")) { dpc.print_e_filt_op = true; dpc.print_simp_filt_op = true; }
        if (vm.count("fzz")) { dpc.run_fzz = true; }

    } catch(std::exception& e) { 

        cout << "Unhandled Exception reached the top of main: " 
                  << e.what() << ", application will now exit" << endl; 
        return false;
    } 

    return true;
}

int main(int argc, char** argv) {

    assert(false);

    /* Process command line */
    DynamicPC dpc;
    std::string posfname(argv[1]);
    Integer max_dim = 2;
    Integer min_time = 0, max_time = std::numeric_limits<Integer>::max();

    if (!parseOptions(argc, argv, posfname, max_dim, min_time, max_time, dpc)) { return -1; }

    std::cout << "posfname: " << posfname << std::endl;
    std::cout << "max dim: " << max_dim << std::endl;
    std::cout << "min/max time: " << min_time << " " << 
        (max_time == std::numeric_limits<Integer>::max() ? std::string("M") : std::to_string(max_time)) 
        << std::endl;


    /* Process dynamic point cloud */

    dpc.init(posfname, max_dim, min_time, max_time);
    std::cout << std::endl << "dpc.init finish" << std::endl;

    dpc.genEvent();
    std::cout << std::endl << "dpc.genEvent finish" << std::endl;

    dpc.initEdgeFilt();
    std::cout << std::endl << "dpc.initEdgeFilt finish" << std::endl;

    dpc.initPersistence();
    std::cout << std::endl << "dpc.initPersistence finish" << std::endl;

    dpc.travEvent();
    std::cout << std::endl << "dpc.travEvent finish" << std::endl;

    std::cout << std::endl << 
        "fw_sw_cnt: " << dpc.fw_sw_cnt << std::endl << 
        "bw_sw_cnt: " << dpc.bw_sw_cnt << std::endl << 
        "ow_sw_cnt: " << dpc.ow_sw_cnt << std::endl << 
        "iw_con_cnt: " << dpc.iw_con_cnt << std::endl << 
        "ow_exp_cnt: " << dpc.ow_exp_cnt << std::endl;

    std::cout << std::endl << "max_e_filt_len: " << dpc.max_e_filt_len << std::endl;
    std::cout << "max_simp_filt_len: " << dpc.max_simp_filt_len << std::endl;

    if (dpc.run_update) 
    { std::cout << std::endl << "total update timing: " << dpc.up_timing.count() << std::endl; }
    if (dpc.run_fzz)
    { std::cout << "total fzz timing: " << dpc.fzz_timing.count() << std::endl; }

    return 0;
}
