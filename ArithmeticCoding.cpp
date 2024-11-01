#include "ArithmeticCoding.hpp"

// Function to build frequency table from input data
void ArithmeticCoding::buildFrequencyTable(const std::string &data, std::map<char, uint32_t> &frequencyTable) {
    for (char ch : data) {
        frequencyTable[ch]++;
    }
}

// Function to build probability table from frequency table
void ArithmeticCoding::buildProbabilityTable(const std::map<char, uint32_t> &frequencyTable, std::map<char, double> &probabilityTable) {
    uint32_t totalSymbols = 0;
    for (const auto &pair : frequencyTable) {
        totalSymbols += pair.second;
    }
    for (const auto &pair : frequencyTable) {
        probabilityTable[pair.first] = static_cast<double>(pair.second) / totalSymbols;
    }
}

// Compresses input file using arithmetic coding
void ArithmeticCoding::compressFile(const std::string &inputFile, const std::string &outputFile) {
    std::ifstream input(inputFile, std::ios::binary);
    std::string data((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    input.close();

    std::map<char, uint32_t> frequencyTable;
    buildFrequencyTable(data, frequencyTable);

    std::map<char, double> probabilityTable;
    buildProbabilityTable(frequencyTable, probabilityTable);

    // Build cumulative probability table
    std::map<char, double> cumulativeProbability;
    double cumulative = 0.0;
    for (const auto &pair : probabilityTable) {
        cumulativeProbability[pair.first] = cumulative;
        cumulative += pair.second;
    }

    // Arithmetic coding
    double low = 0.0;
    double high = 1.0;
    for (char ch : data) {
        double range = high - low;
        high = low + range * (cumulativeProbability[ch] + probabilityTable[ch]);
        low = low + range * cumulativeProbability[ch];
    }

    // Write compressed data to output file
    std::ofstream output(outputFile, std::ios::binary);
    // Write frequency table size
    uint32_t tableSize = frequencyTable.size();
    output.write(reinterpret_cast<char *>(&tableSize), sizeof(tableSize));
    // Write frequency table
    for (const auto &pair : frequencyTable) {
        output.put(pair.first);
        output.write(reinterpret_cast<const char *>(&pair.second), sizeof(pair.second));
    }
    // Write the final low value
    output.write(reinterpret_cast<const char *>(&low), sizeof(low));
    output.close();
}

// Decompresses file compressed with arithmetic coding
void ArithmeticCoding::decompressFile(const std::string &inputFile, const std::string &outputFile) {
    std::ifstream input(inputFile, std::ios::binary);
    if (!input.is_open()) {
        std::cerr << "Failed to open compressed file." << std::endl;
        return;
    }

    // Read frequency table size
    uint32_t tableSize;
    input.read(reinterpret_cast<char *>(&tableSize), sizeof(tableSize));

    // Read frequency table
    std::map<char, uint32_t> frequencyTable;
    for (uint32_t i = 0; i < tableSize; ++i) {
        char ch;
        uint32_t freq;
        input.get(ch);
        input.read(reinterpret_cast<char *>(&freq), sizeof(freq));
        frequencyTable[ch] = freq;
    }

    // Build probability and cumulative probability tables
    std::map<char, double> probabilityTable;
    buildProbabilityTable(frequencyTable, probabilityTable);

    std::map<char, double> cumulativeProbability;
    double cumulative = 0.0;
    for (const auto &pair : probabilityTable) {
        cumulativeProbability[pair.first] = cumulative;
        cumulative += pair.second;
    }

    // Read the final low value
    double codeValue;
    input.read(reinterpret_cast<char *>(&codeValue), sizeof(codeValue));
    input.close();

    // Reconstruct data size
    uint32_t totalSymbols = 0;
    for (const auto &pair : frequencyTable) {
        totalSymbols += pair.second;
    }

    // Arithmetic decoding
    std::string decodedData;
    for (uint32_t i = 0; i < totalSymbols; ++i) {
        double range = 1.0;
        double value = (codeValue - 0.0) / range;

        char symbol = '\0';
        for (const auto &pair : cumulativeProbability) {
            if (value >= pair.second && value < pair.second + probabilityTable[pair.first]) {
                symbol = pair.first;
                break;
            }
        }

        if (symbol == '\0') {
            std::cerr << "Decoding error: symbol not found." << std::endl;
            return;
        }

        decodedData += symbol;

        // Update codeValue for next symbol
        double low = cumulativeProbability[symbol];
        double high = low + probabilityTable[symbol];
        codeValue = (codeValue - low * range) / (high - low);
    }

    // Write decoded data to output file
    std::ofstream output(outputFile, std::ios::binary);
    output.write(decodedData.c_str(), decodedData.size());
    output.close();
}
