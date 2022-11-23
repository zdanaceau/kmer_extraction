#include <iostream>
#include <fstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <chrono>
using namespace std::chrono;


/// @brief Function for extracting kmers from a string
/// @param sequenceBin string containing the sequence to pull kmers from
/// @param ksize size of the kmer
/// @return kmerCollection - map pairing each unique kmer with its frequency
std::map<std::string, int> KmerProfile(std::string sequenceBin, int ksize){
    
    std::map<std::string, int> kmerCollection; 
    for(int i = 0; i < sequenceBin.length() - ksize + 1; i++) {
        
        std::string ithKmer = sequenceBin.substr(i, ksize);

        // Check that the kmer contains valid DNA characters
        bool validKmer = true;
        for(int i = 0; i < ksize; i++){
            // This by default ignores soft masked sequence regions an option to include them can be added
            if(ithKmer[i] != 'A' && ithKmer[i] != 'C'
            && ithKmer[i] != 'T' && ithKmer[i] != 'G') {
                validKmer = false;
            }
        }

        // Count valid kmers appropriately
        if(validKmer){
            auto kmerLoc = kmerCollection.find(ithKmer);

            if(kmerLoc == kmerCollection.end()) {
                kmerCollection.insert(std::pair<std::string,int>(ithKmer, 1));
            }

            else {
                int prevKmerCount = kmerLoc -> second;
                kmerCollection[ithKmer] = prevKmerCount + 1;
            }
        }
    }
    return kmerCollection;
}

/// @brief Function to write kmers from the map specified in KmerProfile to a file
/// @param binLabel Title for sequence bin
/// @param kmerBin map<std::string, int> pairing unique kmer sequences to their frequency
/// @param fileName file to write kmers to
void WriteKmers (std::string binLabel, std::map<std::string, int> kmerBin, std::string fileName) {
    std::ofstream kmerOutput;
    kmerOutput.open(fileName, std::ios::app);
    kmerOutput << binLabel + "\n";
    for(const auto &kmerPair: kmerBin){
        kmerOutput << kmerPair.first << "\t" << std::to_string(kmerPair.second) << "\n";
    }
    kmerOutput.close();
}

/// @brief Main method
/// @param argc 
/// @param argv argv[1] - Input file
///             argv[2] - Kmer size
///             argv[3] - New chromosome marker
///             argv[4] - bin size in bp
/// @return 
int main(int argc, char *argv[])
{

    // print command line arguments
    std::cout << "Input File: " << argv[1] << std::endl;
    std::cout << "Kmer Size: " << argv[2] << std::endl;
    std::cout << "Marker for New Chromosome: " << argv[3] << std::endl;
    std::cout << "Bin size: " << argv[4] << std::endl;

    
    std::string inputFile(argv[1]); // Input file from cmd line argument 1
    int ksize = atoi(argv[2]); // ksize from cmd line argument 2

    // Substring denoting a chromosome 
    std::string newAssemblyMarker(argv[3]); // Chromosome denoting substring from cmd line argument 3
    int binSize = atoi(argv[4]); // Binsize from cmd line argument 4

    std::ifstream inputFasta;
    inputFasta.open(inputFile);
    int nlines = 0;
    // Timer to time operations
    auto totalStart = high_resolution_clock::now();

    // Check file and begin reading
    if(inputFasta.good() && ksize > 0 && binSize > ksize && newAssemblyMarker.length() >= 1) {
        std::string directoryName = inputFile.substr(0, inputFile.find('.')) + "_ksize_" + argv[2] + "_binsize_" + argv[4];
        std::string directoryCommand = "mkdir " + directoryName;
        int returned = system(directoryCommand.c_str());
        if(returned != 0) {
            std::cout << "Output directory already exists, files inside may be overwritten" << std::endl;
        }


        // Initialize parameters for reading sequence
        std::string fastaLine; // Holds current line
        bool readSeq = false; // Dictates whether the seqquence should be read (end of leading 'N')
        bool assemblyStart = false; // Tracks whether an assembly is currently being read
        std::string sequenceBin; // Bin to hold the sequence
        long sequenceIndex = 0; // Current index in the sequence
        long binStart = 0; // Start of the sequence
        long binEnd;
        std::string sequenceID;
        std::ofstream kmerOutput;
        std::string currentWrite;

        while(getline(inputFasta, fastaLine)) {

            // If files contain abnormal newline characters such as \r this may need
            // to be reimplemented in a c++11 compatable format. 
            // This does not appear to be an issue for our assembly files so I'm leaving it out.
            //if(assemblyStart){
                // Erase any possible newline characters to avoid platform dependent conflicts
                //fastaLine.erase(std::remove(fastaLine.begin(), fastaLine.end(), '\n'), fastaLine.cend());
                //fastaLine.erase(std::remove(fastaLine.begin(), fastaLine.end(), '\r'), fastaLine.cend());
            //}


            // Detect when a new assembly sequence starts and starts a new write file
            if(fastaLine.find(newAssemblyMarker)!=std::string::npos && !assemblyStart) {
                //fastaLine.erase(std::remove(fastaLine.begin(), fastaLine.end(), '\n'), fastaLine.cend());
                //fastaLine.erase(std::remove(fastaLine.begin(), fastaLine.end(), '\r'), fastaLine.cend());
                if(fastaLine.length() > 0) {
                    assemblyStart = true;
                    sequenceID = fastaLine;
                }

                // TODO: This should be changed to a better formatted filename
                //currentWrite = directoryName + "/" + "chromosome_" + fastaLine[1] + ".txt";
                // Adjusted to account for larger nums
                currentWrite = directoryName + "/" + "chromosome_" + fastaLine.substr(1, fastaLine.find(' ') - 1) + ".txt";
                kmerOutput.open(currentWrite);
                kmerOutput << "";
                kmerOutput.close();
            }

            // When a new chromosome assembly is detected, end the current assembly and reset to
            // handle the new chromosome sequence
            else if(fastaLine.find(newAssemblyMarker)!=std::string::npos && assemblyStart){
                if(sequenceBin.length() > 0){
                    auto bin = KmerProfile(sequenceBin, ksize);
                    binEnd = binStart + sequenceBin.length() - 1;
                    std::string binLabel = ">" + sequenceID + " Start: " + std::to_string(binStart) + " End: " + std::to_string(binEnd);
                    WriteKmers(binLabel, bin, currentWrite);
                }

                // Reset variables
                binStart = 0;
                readSeq = false;
                sequenceBin = "";
                sequenceIndex = 0;

                // Start writing new assembly
                sequenceID = fastaLine;

                // TODO: This should be changed to a better formatted filename
                //currentWrite = directoryName + "/" + "chromosome_" + fastaLine[1] + ".txt";
                currentWrite = directoryName + "/" + "chromosome_" + fastaLine.substr(1, fastaLine.find(' ') - 1) + ".txt";
                kmerOutput.open(currentWrite);
                kmerOutput << "";
                kmerOutput.close();
            }

            // Handles the case where a new sequence has started but it is not targeted for
            // Analsis
            else if(fastaLine.find(">") != std::string::npos && assemblyStart){
                if(sequenceBin.length() > 0) {
                    auto bin = KmerProfile(sequenceBin, ksize);
                    binEnd = binStart + sequenceBin.length() - 1;
                    std::string binLabel = ">" + sequenceID + " Start: " + std::to_string(binStart) + " End: " + std::to_string(binEnd);
                    WriteKmers(binLabel, bin, currentWrite);
                }
                binStart = 0;
                sequenceBin = "";
                sequenceIndex = 0;
                assemblyStart = false;
            }

            // Start reading sequence and binning after the end of leading N characters
            else if(!readSeq && assemblyStart) {
                int startIndex = sequenceIndex;
                // Determine if the actual read has started yet or if its just another n row
                // for(int index = 0; index < fastaLine.length(); index ++) {
                //     if (fastaLine[index] != 'N'){
                //         readSeq = true;
                //         startIndex = index;
                //         break;
                //     }
                // }

                readSeq = true;
                // Start a 5kb binned string and assign its starting point
                if (readSeq && assemblyStart) {
                    std::string seqToAdd = fastaLine.substr(startIndex);
                    sequenceBin = sequenceBin + fastaLine.substr(startIndex);
                    binStart = sequenceIndex + startIndex;
                }
                sequenceIndex += fastaLine.length();
            }

            // Build each sequenceBin until its full
            else if (sequenceBin.length() + fastaLine.length() < binSize && assemblyStart) {
                sequenceBin = sequenceBin + fastaLine;
                sequenceIndex += fastaLine.length();
            }

            // Construct kmers from bin and write to file 
            // Resets all bin params to begin again 
            else if(assemblyStart){

                // Construct bin and get appropriate lengths 
                int lenToAdd = binSize - sequenceBin.length();
                sequenceBin = sequenceBin + fastaLine.substr(0, lenToAdd);
                int sequenceBinLen = sequenceBin.length();
                binEnd = binStart + sequenceBin.length() - 1;
                auto bin = KmerProfile(sequenceBin, ksize);

                // TODO: Apply compression function if desired
                // Create kmer label and write to file
                std::string binLabel = ">" + sequenceID + " Start: " + std::to_string(binStart) + " End: " + std::to_string(binEnd);
                WriteKmers(binLabel, bin, currentWrite);
                
                // Reset bin params 
                binStart = binEnd + 1;
                sequenceBin = fastaLine.substr(lenToAdd, fastaLine.length() - lenToAdd);
                int testBinLength = sequenceBin.length();
                sequenceIndex += fastaLine.length();

            }

        }
        // Construct final bin 
        if(sequenceBin.length() > 0){
            auto bin = KmerProfile(sequenceBin, ksize);
            binEnd = binStart + sequenceBin.length();
            std::string binLabel = ">" + sequenceID + " Start: " + std::to_string(binStart) + " End: " + std::to_string(binEnd);
            WriteKmers(binLabel, bin, currentWrite);
        }
    }
    else{
        std::cout<<"Invalid input file or ksize"<<std::endl;
    }
    inputFasta.close();
    auto totalStop = high_resolution_clock::now();
    auto totalDuration = duration_cast<microseconds>(totalStop - totalStart);
    std::cout << "Duration us: " << totalDuration.count() << std::endl;
}