#ifndef ARITHMETICCODING_HPP
#define ARITHMETICCODING_HPP

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cstdint>
#include <string>

class ArithmeticCoding {
public:
    void compressFile(const std::string &inputFile, const std::string &outputFile);
    void decompressFile(const std::string &inputFile, const std::string &outputFile);

private:
    void buildFrequencyTable(const std::string &data, std::map<char, uint32_t> &frequencyTable);
    void buildProbabilityTable(const std::map<char, uint32_t> &frequencyTable, std::map<char, double> &probabilityTable);
};

#endif // ARITHMETICCODING_H
