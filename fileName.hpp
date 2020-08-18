# pragma once
# include <iomanip>
# include <string>

// Change to_string with precision for main folder name
template <typename T>
std::string to_stringWithPrecision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::fixed <<std::setprecision(n) << a_value;
    return out.str();
}

template <typename T>
std::string to_stringWithExponent(const T &a_value, const int &n=1)
{
    std::ostringstream out;
    out<<std::scientific<<std::setprecision(n) << a_value;
    return out.str();
}


std::string fileName(const int &t_networkSize, const double &t_alpha, const double &t_gamma, const int &t_ensembleSize){
    std::string fileName="N="+to_stringWithExponent((double)t_networkSize,1)+",a="+to_stringWithPrecision(t_alpha,1)+",r="+to_stringWithPrecision(t_gamma,1)+",e="+std::to_string(t_ensembleSize);
    return fileName;
}

std::string fileName(const int &t_networkSize, const double &t_g, const int &t_ensembleSize){
    std::string fileName="N="+to_stringWithExponent((double)t_networkSize,1)+",g="+to_stringWithPrecision(t_g,1)+",e="+std::to_string(t_ensembleSize);
    return fileName;
}

std::string fileName(const int &t_networkSize, const double &t_g){
    std::string fileName="N="+to_stringWithExponent((double)t_networkSize,1)+",g="+to_stringWithPrecision(t_g,1);
    return fileName;
}