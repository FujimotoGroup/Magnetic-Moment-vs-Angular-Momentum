#include <iostream>
#include <iomanip>
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/optional.hpp>


void load_ini() {
    using namespace boost::property_tree;

    ptree pt;
    read_ini("./ini/parameter.ini", pt);

    // physics parameters {{{
    if (boost::optional<double> a = pt.get_optional<double>("physics.angstrom")) {
        angstrom = a.get();
    } else {
        std::cout << "'hbar' is nothing" << std::endl;
    }

    if (boost::optional<double> h = pt.get_optional<double>("physics.hbar")) {
        hbar = h.get();
    } else {
        std::cout << "'hbar' is nothing" << std::endl;
    }

    if (boost::optional<double> m = pt.get_optional<double>("physics.mass")) {
        mass = m.get();
    } else {
        std::cout << "'angstrom' is nothing" << std::endl;
    }

    if (boost::optional<double> c = pt.get_optional<double>("physics.charge")) {
        charge = c.get();
    } else {
        std::cout << "'hbar' is nothing" << std::endl;
    }
    // }}}

    // settings {{{
    if (boost::optional<int> b = pt.get_optional<int>("settings.bands")) {
        bands = b.get();
    } else {
        std::cout << "'bands' is nothing" << std::endl;
    }

    if (boost::optional<int> bT = pt.get_optional<int>("settings.bandsT")) {
        bandsT = bT.get();
    } else {
        std::cout << "'bandsT' is nothing" << std::endl;
    }

    if (boost::optional<int> bL = pt.get_optional<int>("settings.bandsL")) {
        bandsL = bL.get();
    } else {
        std::cout << "'bandsL' is nothing" << std::endl;
    }
    // }}}
}

