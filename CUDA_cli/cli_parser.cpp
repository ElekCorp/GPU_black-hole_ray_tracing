#include "cli_parser.h"

void parse_args(int argc, char* argv[], Params& p)
{
    std::unordered_map<std::string, double*> float_args = {
        {"--errormax", &p.errormax},
        {"--de0", &p.de0},
        {"--rs", &p.rs},
        {"--delta_a", &p.delta_a},
        {"--a", &p.a},
        {"--Q", &p.Q},
        {"--t0", &p.t_0},
        {"--r0", &p.r_0},
        {"--theta0", &p.theta_0},
        {"--phi0", &p.phi_0},
        {"--kepernyo_high", &p.kepernyo_high},
        {"--kepernyo_tav", &p.kepernyo_tav},
        {"--sugar_ki", &p.sugar_ki},
        {"--gyuru_kicsi", &p.gyuru_sugar_kicsi},
        {"--gyuru_nagy", &p.gyuru_sugar_nagy},
    };

    std::unordered_map<std::string, int*> int_args = {
        {"--kepernyoSZELES", &p.kepernyoSZELES},
        {"--kepernyoMAGAS", &p.kepernyoMAGAS},
        {"--SZELES", &p.SZELES},
        {"--MAGAS", &p.MAGAS},
        {"--ikezd", &p.ikezd},
        {"--jkezd", &p.jkezd},
        {"--iveg", &p.iveg},
    };

    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (float_args.count(arg)) {
            *float_args[arg] = std::stod(argv[++i]);
        }
        else if (int_args.count(arg)) {
            *int_args[arg] = std::stoi(argv[++i]);
        }
    }
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "--double") p.prec = Precession::Double;
        if (arg == "--float")  p.prec = Precession::Float;
    }

    if(p.iveg == 0)
    {
        p.iveg=p.kepernyoSZELES;
    }
}

