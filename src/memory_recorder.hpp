#ifndef memory_recorder_h_
#define memory_recorder_h_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <ctime>
#include <map>
#include <sys/stat.h>
#include <numeric>
#include <mpi.h>
#include <sys/resource.h>
#include <unistd.h>

struct MemoryRecorder
{
    const long kb = 1024;
    const long mb = 1024*1024;

    int getrss = 0, getrss_summary = 0, getmeminfo = 0;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    int globalrank, localrank, localsize;
    int pid;

    MPI_Comm shmcomm;

    long *rss_collect;
    long maxrss, local_maxrss, global_maxrss;

    std::string dirOut, memfileOut, rssfileOut;
    std::map<std::string, std::vector<long>> freemem;
    std::vector<std::string> meminfo_names;

    MemoryRecorder();
    void getMaxRSS();
    void read_meminfo(std::string const &loc);
    void summarizeMaxRSS();
    void write_meminfo();
    void write_rss();
    ~MemoryRecorder();
};

MemoryRecorder::MemoryRecorder() {

    int k, namelen;
    MPI_Comm_rank(MPI_COMM_WORLD, &globalrank);
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0,
        MPI_INFO_NULL, &shmcomm);
    MPI_Comm_rank(shmcomm, &localrank);
    MPI_Comm_size(shmcomm, &localsize);
    MPI_Get_processor_name(hostname, &namelen);
    pid = getpid();
    std::string strhost;
    for (k=0; k<namelen; k++) {
        strhost += hostname[k];
    }

    // Get Memory Files to read from.
    k=0;
    char fname_in[128];
    std::ifstream stream;
    stream.open("/usr/bin/bash", std::ios::in);
    while (stream.good()) {
        stream.close();
        snprintf(fname_in, sizeof(fname_in),
            "/sys/devices/system/node/node%d/meminfo",k);
        stream.open(fname_in, std::ios::in);
        k++;
        meminfo_names.push_back(fname_in);
    }
    meminfo_names.pop_back();

    // Set paths
    dirOut = "memory_records";
    memfileOut = "meminfo_" + strhost + ".csv";
    rssfileOut = "rss_" + strhost + ".out";
    rss_collect = (long*)malloc(sizeof(long) * localsize);
}

void MemoryRecorder::getMaxRSS() {
    getrss++; 
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    maxrss = usage.ru_maxrss; // In kb
}

void MemoryRecorder::summarizeMaxRSS() {
    getrss_summary = 1;
    if ( getrss == 0 ) {
        this->getMaxRSS();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Reduce(&maxrss, &global_maxrss, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&maxrss, &local_maxrss, 1, MPI_LONG, MPI_SUM, 0, shmcomm);
    MPI_Gather(&maxrss, 1, MPI_LONG, rss_collect, 1, MPI_LONG, 0, shmcomm);
}

void MemoryRecorder::read_meminfo(std::string const &loc) {
    getmeminfo++;
    if (localrank == 0) {
        std::vector<long> memnow;
        std::string line;
        std::string word;
        std::ifstream stream;

        for (std::string fn : meminfo_names) {
            stream.open(fn, std::ios::in);
            // stream.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
            std::getline(stream,line);
            // Extract the fourth field from the second line
            for (int i = 0; i < 4; i++) {
                stream >> word; // Reads one space separated word from the stream.
                // std::cout << word << " ";
            }
            memnow.push_back(std::stol(word));
        }

        long memsum = std::accumulate(memnow.begin(), memnow.end(), 0);
        memnow.push_back(memsum);
        freemem.insert({loc, memnow});
    }
    return;
}

void MemoryRecorder::write_meminfo() {
    
    if ( getmeminfo == 0 ) {
        if (globalrank == 0) {
            std::cout << "No Meminfo readings to write out " << std::endl;
        }
        return;
    }

    if (localrank == 0) {

        // Create the directory if it does not exist
        struct stat st = {0};
        if (stat(dirOut.c_str(), &st) == -1) {
            mkdir(dirOut.c_str(), 0755);
        }
        // Open the file in the directory
        std::ofstream outfile;
        outfile.open(dirOut + "/" + memfileOut, std::ios::out);
        outfile << "code_location,";
        for (int i=0; i<meminfo_names.size(); i++) {
            if (i != meminfo_names.size()-1) {
                outfile << i << ",";
            } else {
                outfile << i << std::endl;
            }
        }
        for (const auto& item : freemem){
            outfile << item.first << ",";
            for ( auto freem : item.second ) {
                if (freem != item.second.back()) {
                    outfile << freem/mb << ",";
                } else {
                    outfile << freem/mb << std::endl;
                }
            }
        }

        outfile.close();
    }

}

void MemoryRecorder::write_rss() {
    if (getrss_summary == 0) {
        this->summarizeMaxRSS();
    }
    if (localrank == 0) {
        struct stat st = {0};
        if (stat(dirOut.c_str(), &st) == -1) {
            mkdir(dirOut.c_str(), 0755);
        }
        // Open the file in the directory
        std::ofstream outfile;
        outfile.open(dirOut + "/" + rssfileOut, std::ios::out);
        outfile << "Node Rss Sum: " << local_maxrss/mb << " (GiB)" << std::endl;
        for (int i=0; i<localsize; i++) {
            outfile << "Rank: " << i << " MaxRSS: " << rss_collect[i]/kb << " (MiB)" << std::endl;
        }
        outfile.close();
    }
    if (globalrank == 0) {
        std::cout << "TOTAL RSS MAX SUM: " << global_maxrss/mb << " (GiB)" << std::endl;
    }
}


MemoryRecorder::~MemoryRecorder() {
    free(rss_collect);
}

#endif


