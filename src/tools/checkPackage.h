#include <cstring>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

class PackageChecker
{

public:
  static string PIP_FREEZE;
  static string PYTHON_VERSION;
  static string R_VERSION;
  static string R_PACKAGE;
  static string PYTHON;
  static string R;
  static string R_PACKAGES;
  PackageChecker();

    struct Package
    {
        std::string name;    // Package name
        std::string version; // Version number. String because it could be something like 12.0.1.
    };

    static string exec(const char *cmd, bool is_installation);

    static string trim(const string &str);

    static bool check_version(const string pkg_version, const string installed_version);

    static int *get_major_minor(const string version);

    static bool doSegment(char *sentence, vector<Package> rpackages, string package_delim);

    static bool check_packages(const std::vector<Package> &, const std::string &, const std::string &);

    static bool check_installation(const std::string &);
    
    static bool checkRAnaquin(const std::string &v)
    {
        std::vector<PackageChecker::Package> p;
        p.push_back({"Anaquin", v});
        return PackageChecker::check_packages(p, PackageChecker::R_PACKAGES, ",");
    }
    
    static bool checkR()
    {
        return PackageChecker::check_installation(PackageChecker::R);
    }

    static bool checkPython()
    {
        return PackageChecker::check_installation(PackageChecker::PYTHON);
    }
};
