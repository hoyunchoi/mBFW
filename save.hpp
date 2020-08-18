#pragma once

#include <filesystem>

#include "/pds/pds172/hoyun1009/library/CSV.hpp"
#include "parameters.hpp"
#include "fileName.hpp"

const std::string folder="/pds/pds172/hoyun1009/mBFW/data/";

void saveeT_c_minus(const std::string &t_filename, const int &t_coreNum, const std::vector<double> &t_eT_c_minus){
    const std::string file=folder+"T_a/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    writeCSV(file, t_eT_c_minus);
}

void saveOrderParameter(const std::string &t_filename, const int &t_coreNum, const std::vector<double> &t_orderParameter){
    const std::string file=folder+"orderParameter/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    writeCSV(file,t_orderParameter);
}

void saveSecondGiant(const std::string &t_filename, const int &t_coreNum, const std::vector<double> &t_secondGiant){
    const std::string file=folder+"secondGiant/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    writeCSV(file,t_secondGiant);
}

void saveOrderParameterAfter(const std::string &t_filename, const int &t_networkSize, const int &t_ensembleSize, const int &t_coreNum, const std::vector<double> &t_orderParameterAfter, const std::vector<int> &t_sampledOrderParameterAfter){
    const std::string file=folder+"orderParameterAfter/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    const std::string file2=folder+"orderParameterAfter/sampled/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    writeCSV(file,t_orderParameterAfter);
    writeCSV(file2, t_sampledOrderParameterAfter);
}

void saveMeanClusterSize(const std::string &t_filename, const int &t_coreNum, const std::vector<double> &t_meanClusterSize){
    const std::string file=folder+"meanClusterSize/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    writeCSV(file,t_meanClusterSize);
}


void saveClusterSizeDistribution(const std::string &t_filename, const int &t_coreNum, const std::vector<double> &t_orderParameterOfClusterSizeDistribution, const std::vector<std::vector<double>> &t_clusterSizeDistribution, const std::vector<int> &t_sampledClusterSizeDistribution, const int &t_exclude){
    const int N=t_clusterSizeDistribution[0].size();
    std::map<double, int> sampled;
    for (int i=0; i<t_orderParameterOfClusterSizeDistribution.size(); ++i){
        sampled[t_orderParameterOfClusterSizeDistribution[i]]=t_sampledClusterSizeDistribution[i];
        const std::string file=folder+"clusterSizeDistribution"+std::to_string(t_exclude)+"/"+t_filename+",m="+to_stringWithPrecision(t_orderParameterOfClusterSizeDistribution[i],4)+"-"+std::to_string(t_coreNum)+".txt";
        std::map<int, double> clusterSizeDistribution;
        for (int j=0; j<N; ++j){
            if (t_clusterSizeDistribution[i][j]!=0){
                clusterSizeDistribution[j]=t_clusterSizeDistribution[i][j];
            }
        }
        writeCSV(file, clusterSizeDistribution);
    }
    const std::string file=folder+"clusterSizeDistribution"+std::to_string(t_exclude)+"/sampled/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    writeCSV(file, sampled);
}

void saveInterEventTimeDistribution(const std::string &t_filename, const int &t_coreNum, const std::vector<std::vector<double>> &t_interEventTimeDistribution){
    const std::vector<std::string> type={"before","after"};
    for (int i=0; i<2; ++i){
        const std::string file=folder+"interEventTimeDistribution/"+type[i]+"/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
        std::map<int, double> interEventTime;
        for (int j=0; j<t_interEventTimeDistribution[i].size(); ++j){
            if (t_interEventTimeDistribution[i][j]!=0){
                interEventTime[j]=t_interEventTimeDistribution[i][j];
            }
        }
        writeCSV(file, interEventTime);
    }
}

void saveOrderParameterDistribution(const std::string &t_filename, const int &t_coreNum, const std::vector<double> &t_timeOfOrderOarameterDistribution, const std::vector<std::vector<double>> &t_orderParameterDistribution){
    const int N=t_orderParameterDistribution[0].size();
    for (int i=0; i<t_timeOfOrderOarameterDistribution.size(); ++i){
        const std::string file=folder+"orderParameterDistribution/"+t_filename+",t="+to_stringWithPrecision(t_timeOfOrderOarameterDistribution[i],4)+"-"+std::to_string(t_coreNum)+".txt";
        std::map<double, double> orderParameterDistribution;
        for (int j=0; j<N; ++j){
            if (t_orderParameterDistribution[i][j]!=0){
                orderParameterDistribution[j/(double)N]=t_orderParameterDistribution[i][j];
            }
        }
        writeCSV(file,orderParameterDistribution);
    }
}

void saveTime(const std::vector<std::vector<double>> &t_mk, const std::vector<std::vector<double>> &t_tu){
    const std::string file1=folder+"mk.txt";
    const std::string file2=folder+"tu.txt";
    writeCSV(file1, t_mk);
    writeCSV(file2, t_tu);
}

void saveLUMK(const std::vector<std::vector<int>> &t_LUMK){
    const std::string file=folder+"LU_PLPU_MK.txt";
    writeCSV(file, t_LUMK);
}

void saveLUPLPU(const std::string &t_filename, const std::vector<std::vector<int>> &t_LUPLPU){
    const std::string file=folder+"LUPLPU/"+t_filename+".txt";
    writeCSV(file, t_LUPLPU);
}

void saveInterEventTime(const std::string &t_filename, const int &t_coreNum, const std::vector<double> &t_interEventTime){
    const std::string file=folder+"interEventTime/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    std::map<int,double> interEventTime;
    for (int i=0; i<t_interEventTime.size(); ++i){
        if (t_interEventTime[i]!=0){
            interEventTime[i]=t_interEventTime[i];
        }
    }
    writeCSV(file, interEventTime);
}

void saveDeltaM(const std::string &t_filename, const int &t_coreNum, const std::vector<std::vector<double>> &t_deltaMDistribution){
    const std::vector<std::string> type={"before","after"};
    for (int i=0; i<2; ++i){
        const std::string file=folder+"deltaMDistribution/"+type[i]+"/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
        std::map<int, double> deltaM;
        for (int j=0; j<t_deltaMDistribution[i].size(); ++j){
            if (t_deltaMDistribution[i][j]!=0){
                deltaM[j]=t_deltaMDistribution[i][j];
            }
        }
        writeCSV(file, deltaM);
    }
}

void saveAcceptanceDistribution(const std::string &t_filename, const int &t_coreNum, const std::vector<std::vector<double>> &t_acceptance){
    const std::string file=folder+"acceptanceDistribution/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    const int precision=10;
    writeCSV(file, t_acceptance, precision);
}

void saveAgeDistribution(const std::string &t_filename, const int &t_coreNum, const std::vector<std::vector<double>> &t_ageDistribution){
    const std::vector<std::string> type={"before","after"};
    for (int i=0; i<2; ++i){
        const std::string file=folder+"ageDistribution/"+type[i]+"/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
        std::map<int, double> age;
        for (int j=1; j<t_ageDistribution[i].size(); ++j){
            if (t_ageDistribution[i][j]!=0){
                age[j]=t_ageDistribution[i][j];
            }
        }
        writeCSV(file, age);
    }
}

void saveK_DeltaAcceptance(const std::string &t_filename, const int &t_coreNum, const std::vector<double> &t_k_DeltaAcceptance){
    const std::string file=folder+"k_DeltaAcceptance/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    std::map<int, double> deltaAcceptance;
    for (int i=0; i<t_k_DeltaAcceptance.size(); ++i){
        if (t_k_DeltaAcceptance[i]!=0){
            deltaAcceptance[i]=t_k_DeltaAcceptance[i];
        }
    }
    const int precision=10;
    writeCSV(file, deltaAcceptance, precision);
}

void saveDeltaK_DeltaAcceptance(const std::string &t_filename, const int &t_coreNum, const std::vector<double> &t_deltaK_DeltaAcceptance){
    const std::string file=folder+"deltaK_DeltaAcceptance/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    std::map<int, double> deltaAcceptance;
    for (int i=0; i<t_deltaK_DeltaAcceptance.size(); ++i){
        if (t_deltaK_DeltaAcceptance[i]!=0){
            deltaAcceptance[i]=t_deltaK_DeltaAcceptance[i];
        }
    }
    const int precision=10;
    writeCSV(file, deltaAcceptance, precision);
}

void saveTime_DeltaAcceptance(const std::string &t_filename, const int &t_coreNum, const std::vector<double> &t_time_DeltaAcceptance){
    const std::string file=folder+"time_DeltaAcceptance/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    std::map<int, double> deltaAcceptance;
    for (int i=0; i<t_time_DeltaAcceptance.size(); ++i){
        if (t_time_DeltaAcceptance[i]!=0){
            deltaAcceptance[i]=t_time_DeltaAcceptance[i];
        }
    }
    const int precision=10;
    writeCSV(file, deltaAcceptance, precision);
}


void savePeriodAcceptanceAreaDistribution(const std::string &t_filename, const int &t_coreNum, const std::vector<std::vector<double>> &t_periodAcceptanceAreaDistribution){
    const std::string file=folder+"periodAcceptanceAreaDistribution/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    const int precision=10;
    writeCSV(file, t_periodAcceptanceAreaDistribution, precision);
}


void saveSuccessiveK(const std::vector<std::vector<int>> &t_successiveK){
    const std::string file=folder+"ssuccessiveK.txt";
    writeCSV(file,t_successiveK,10);
}

void savePeriodAcceptance_UpperBoundRatio(const std::string &t_filename, const int &t_coreNum, const std::vector<std::vector<double>> &t_periodAcceptance_UpperBoundRatio){
    const std::string file=folder+"periodAcceptance_UpperBoundRatio/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    const int precision=10;
    writeCSV(file, t_periodAcceptance_UpperBoundRatio, precision);
}


void savePeriodAcceptance_DeltaK(const std::string &t_filename, const int &t_coreNum, const std::vector<std::vector<double>> &t_periodAcceptance_DeltaK){
    const std::string file=folder+"periodAcceptance_DeltaK/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    const int precision=10;
    writeCSV(file, t_periodAcceptance_DeltaK, precision);
}


void saveInterEventTimeAcceptance(const std::string &t_filename, const int &t_coreNum, const std::vector<double> &t_interEventTime_Acceptance){
    const std::string file=folder+"interEventTime_Acceptance/"+t_filename+"-"+std::to_string(t_coreNum)+".txt";
    std::map<int, double> interEventTime_Acceptance;
    for (int i=0; i<t_interEventTime_Acceptance.size(); ++i){
        if(t_interEventTime_Acceptance[i]){
            interEventTime_Acceptance[i]=t_interEventTime_Acceptance[i];
        }
    }
    const int precision=10;
    writeCSV(file, interEventTime_Acceptance, precision);
}