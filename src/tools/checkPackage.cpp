#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <assert.h>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>
#include <algorithm>
#include <regex>
#include "checkPackage.h"

string PackageChecker::PIP_FREEZE = "pip freeze";
string PackageChecker::PYTHON_VERSION = "python --version";
string PackageChecker::R_VERSION = "R --version";
string PackageChecker::R_PACKAGES = "R --salve -e 'a<-as.data.frame(installed.packages()); a$pac_ver <- paste(a$Package, a$Version, sep=\",\"); print(as.data.frame(a$pac_ver), row.names=FALSE, column.names=FALSE)'";

string PackageChecker::PYTHON = "which python";
string PackageChecker::R = "which R";

string PackageChecker::exec(const char *cmd, bool is_installation)
{
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe)
    {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()))
    {
        // R package regex
        std::regex regexR("[a-zA-Z0-9.-]+,[0-9.-]+");
        // Python package regex
        std::regex regexPy("[a-zA-Z0-9.-]+==[0-9.-]+");
        std::string buffer_data = PackageChecker::trim(buffer.data());
        if (is_installation || std::regex_search(buffer_data, regexR) || std::regex_search(buffer_data, regexPy))
        {
            result += buffer_data;
        }
    }
    return result;
}

string PackageChecker::trim(const string &str)
{
    size_t first = str.find_first_not_of(' ');
    if (string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

bool PackageChecker::check_version(const string pkg_version, const string installed_version)
{
    int *pkg_major_minor = PackageChecker::get_major_minor(pkg_version);
    int *installed_major_minor = PackageChecker::get_major_minor(installed_version);

    // cout << "Package_version" << pkg_major_minor[0] << pkg_major_minor[1] << endl;
    // cout << "install version" << installed_major_minor[0] << installed_major_minor[1] << endl;

    if (installed_major_minor[0] >= pkg_major_minor[0])
    {
        if (installed_major_minor[1] >= pkg_major_minor[1])
        {
            return true;
        }
    }

    return false;
}

int *PackageChecker::get_major_minor(const string version)
{
    std::regex rx("([0-9]+)\\.([0-9]+)");
    std::match_results<std::string::const_iterator> mr;

    std::regex_search(version, mr, rx);

    int *major_minor = new int[2];

    if (mr.size() > 2)
    {
        major_minor[0] = stoi(mr[1]); // major
        major_minor[1] = stoi(mr[2]); // minor
    }
    else
    {
        major_minor[0] = 0; // major
        major_minor[1] = 0; // minor
    }
    return major_minor;
}

bool PackageChecker::doSegment(char *sentence, vector<Package> rpackages, string package_delim)
{
    for (auto x : rpackages)
    {
        std::stringstream ss(sentence);
        std::string to;
        std::string wo;
        bool break_first_while = false;

        while (std::getline(ss, to, '\n'))
        {
            to = PackageChecker::trim(to);
            if (to.find(package_delim) != std::string::npos)
            {
                string installed_package = to.substr(0, to.find(package_delim));
                string installed_version = to.substr(to.find(package_delim) + package_delim.length(), to.length());

                if (x.name == installed_package)
                {
                    if (PackageChecker::check_version(x.version, installed_version))
                    {
                        break_first_while = true;
                        break;
                    }
                    else
                    {
                        return false;
                    }
                }
            }
            else
            {
                continue;
            }

            if (break_first_while)
            {
                break;
            }
        }

        if (!break_first_while)
        {
            return false;
        }
    }

    return true;
}

bool PackageChecker::check_packages(const std::vector<Package> &packages, const std::string &packages_cmd, const std::string &package_delim)
{
    const auto versions = PackageChecker::exec(packages_cmd.c_str(), false);
    int n = versions.length();
    char char_array[n + 1];
    return PackageChecker::doSegment(strcpy(char_array, versions.c_str()), packages, package_delim);
}

bool PackageChecker::check_installation(const std::string &name)
{
    return !PackageChecker::exec(name.c_str(), true).empty();
}
