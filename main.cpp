#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <queue>
#include <stdexcept>p
#include <stdint.h>

// Helper to write bits to file
class BitWriter {
    std::ofstream &out;
    uint8_t buffer;
    int bit_count;

public:
    explicit BitWriter(std::ofstream &outFile) : out(outFile), buffer(0), bit_count(0) {}

    void writeBit(bool bit) {
        buffer |= bit << (7 - bit_count);
        ++bit_count;
        if (bit_count == 8) {
            out.put(buffer);
            bit_count = 0;
            buffer = 0;
        }
    }

    void writeBits(uint32_t bits, int length) {
        for (int i = length - 1; i >= 0; --i) {
            writeBit((bits >> i) & 1);
        }
    }

    void flush() {
        if (bit_count > 0) {
            out.put(buffer);
            bit_count = 0;
            buffer = 0;
        }
    }
};

// Helper to read bits from file
class BitReader {
    std::ifstream &in;
    uint8_t buffer;
    int bit_count;

public:
    explicit BitReader(std::ifstream &inFile) : in(inFile), buffer(0), bit_count(0) {}

    bool readBit() {
        if (bit_count == 0) {
            buffer = in.get();
            if (in.eof()) {
                throw std::runtime_error("Unexpected EOF in BitReader");
            }
            bit_count = 8;
        }
        bool bit = (buffer >> (bit_count - 1)) & 1;
        --bit_count;
        return bit;
    }

    uint32_t readBits(int length) {
        uint32_t value = 0;
        for (int i = 0; i < length; ++i) {
            value = (value << 1) | readBit();
        }
        return value;
    }
};

// Efficient Burrows-Wheeler Transform using Suffix Array
std::pair<std::string, int> burrows_wheeler_transform(const std::string &input) {
    int n = input.size();
    std::vector<int> suffix_array(n);

    // Initialize suffix array
    for (int i = 0; i < n; ++i) {
        suffix_array[i] = i;
    }

    // Sort suffixes
    std::sort(suffix_array.begin(), suffix_array.end(), [&input, n](int a, int b) {
        return strcmp(input.c_str() + a, input.c_str() + b) < 0;
    });

    // Build BWT result from last characters of sorted suffixes
    std::string bwt;
    bwt.reserve(n);
    int primary_index = -1;
    for (int i = 0; i < n; ++i) {
        if (suffix_array[i] == 0) {
            primary_index = i;
            bwt += input[n - 1];
        } else {
            bwt += input[suffix_array[i] - 1];
        }
    }

    return {bwt, primary_index};
}

// Inverse Burrows-Wheeler Transform
std::string inverse_burrows_wheeler_transform(const std::string &bwt, int primary_index) {
    int n = bwt.size();
    std::vector<int> count(256, 0);
    std::vector<int> next(n);

    // Count occurrences
    for (char c : bwt) {
        ++count[(unsigned char)c];
    }

    // Cumulative counts
    std::vector<int> cumulative_count(256, 0);
    int sum = 0;
    for (int i = 0; i < 256; ++i) {
        cumulative_count[i] = sum;
        sum += count[i];
    }

    // Build next array
    for (int i = 0; i < n; ++i) {
        char c = bwt[i];
        next[cumulative_count[(unsigned char)c]++] = i;
    }

    // Reconstruct the original string
    std::string original;
    original.resize(n);
    int idx = next[primary_index];
    for (int i = n - 1; i >= 0; --i) {
        original[i] = bwt[idx];
        idx = next[idx];
    }

    return original;
}

// Move-to-Front Encoding using linked list for efficiency
std::vector<int> move_to_front_encode(const std::string &bwt) {
    std::vector<int> mtfEncoded;
    mtfEncoded.reserve(bwt.size());

    std::vector<uint8_t> symbolTable(256);
    for (int i = 0; i < 256; ++i) {
        symbolTable[i] = i;
    }

    for (unsigned char c : bwt) {
        int index = 0;
        while (symbolTable[index] != c) {
            ++index;
        }
        mtfEncoded.push_back(index);

        // Move to front
        while (index > 0) {
            symbolTable[index] = symbolTable[index - 1];
            --index;
        }
        symbolTable[0] = c;
    }

    return mtfEncoded;
}

// Move-to-Front Decoding
std::string move_to_front_decode(const std::vector<int> &mtf) {
    std::string decoded;
    decoded.reserve(mtf.size());

    std::vector<uint8_t> symbolTable(256);
    for (int i = 0; i < 256; ++i) {
        symbolTable[i] = i;
    }

    for (int index : mtf) {
        unsigned char c = symbolTable[index];
        decoded += c;

        // Move to front
        while (index > 0) {
            symbolTable[index] = symbolTable[index - 1];
            --index;
        }
        symbolTable[0] = c;
    }

    return decoded;
}

// Huffman Tree node
struct HuffmanNode {
    int symbol;
    uint64_t freq;
    HuffmanNode *left;
    HuffmanNode *right;

    HuffmanNode(int s, uint64_t f) : symbol(s), freq(f), left(nullptr), right(nullptr) {}
};

// Comparator for the priority queue
struct Compare {
    bool operator()(HuffmanNode *left, HuffmanNode *right) {
        return left->freq > right->freq;
    }
};

// Huffman Encoding
class HuffmanEncoder {
    std::map<int, std::vector<bool>> huffmanCode;

    void buildCode(HuffmanNode *root, std::vector<bool> &code) {
        if (!root) return;

        if (root->symbol != -1) {
            huffmanCode[root->symbol] = code;
        } else {
            code.push_back(false);
            buildCode(root->left, code);
            code.pop_back();

            code.push_back(true);
            buildCode(root->right, code);
            code.pop_back();
        }
    }

public:
    explicit HuffmanEncoder(const std::map<int, uint64_t> &freq) {
        // Build Huffman Tree
        std::priority_queue<HuffmanNode *, std::vector<HuffmanNode *>, Compare> pq;

        for (auto pair : freq) {
            pq.push(new HuffmanNode(pair.first, pair.second));
        }

        if (pq.empty()) {
            return;
        }

        while (pq.size() > 1) {
            HuffmanNode *left = pq.top();
            pq.pop();
            HuffmanNode *right = pq.top();
            pq.pop();

            auto *newNode = new HuffmanNode(-1, left->freq + right->freq);
            newNode->left = left;
            newNode->right = right;
            pq.push(newNode);
        }

        std::vector<bool> code;
        buildCode(pq.top(), code);
    }

    const std::map<int, std::vector<bool>> &getHuffmanCode() const {
        return huffmanCode;
    }
};

// Huffman Decoding
class HuffmanDecoder {
    HuffmanNode *root;

public:
    explicit HuffmanDecoder(const std::map<int, uint64_t> &freq) {
        // Build Huffman Tree
        std::priority_queue<HuffmanNode *, std::vector<HuffmanNode *>, Compare> pq;

        for (auto pair : freq) {
            pq.push(new HuffmanNode(pair.first, pair.second));
        }

        if (pq.empty()) {
            root = nullptr;
            return;
        }

        while (pq.size() > 1) {
            HuffmanNode *left = pq.top();
            pq.pop();
            HuffmanNode *right = pq.top();
            pq.pop();

            auto *newNode = new HuffmanNode(-1, left->freq + right->freq);
            newNode->left = left;
            newNode->right = right;
            pq.push(newNode);
        }

        root = pq.top();
    }

    std::vector<int> decode(BitReader &reader, uint64_t numSymbols) {
        std::vector<int> decoded;
        decoded.reserve(numSymbols);

        if (!root) {
            return decoded;
        }

        for (uint64_t i = 0; i < numSymbols; ++i) {
            HuffmanNode *node = root;
            while (node->left || node->right) {
                bool bit;
                try {
                    bit = reader.readBit();
                } catch (const std::runtime_error &e) {
                    throw std::runtime_error("Error during Huffman decoding: " + std::string(e.what()));
                }
                node = bit ? node->right : node->left;
                if (!node) {
                    throw std::runtime_error("Invalid Huffman code encountered during decoding.");
                }
            }
            decoded.push_back(node->symbol);
        }
        return decoded;
    }
};

// Compression Function
void compress(const std::string &inputFile, const std::string &outputFile) {
    // Step 1: Read input file
    std::ifstream inFile(inputFile, std::ios::binary);
    if (!inFile) {
        std::cerr << "Cannot open input file: " << inputFile << std::endl;
        return;
    }
    std::string input((std::istreambuf_iterator<char>(inFile)), std::istreambuf_iterator<char>());
    inFile.close();

    // Handle empty input
    if (input.empty()) {
        std::ofstream outFile(outputFile, std::ios::binary);
        uint64_t numSymbols = 0;
        outFile.write(reinterpret_cast<const char *>(&numSymbols), sizeof(numSymbols));
        for (int i = 0; i < 256; ++i) {
            uint64_t frequency = 0;
            outFile.write(reinterpret_cast<const char *>(&frequency), sizeof(frequency));
        }
        outFile.close();
        return;
    }

    // Step 2: Efficient Burrows-Wheeler Transform
    auto [bwt, primary_index] = burrows_wheeler_transform(input);

    // Step 3: Efficient Move-to-Front Encoding
    std::vector<int> mtfEncoded = move_to_front_encode(bwt);

    // Step 4: Compute frequencies
    std::map<int, uint64_t> freq;
    for (int symbol : mtfEncoded) {
        freq[symbol]++;
    }

    // Step 5: Build Huffman tree
    HuffmanEncoder huffman(freq);
    const auto &huffmanCode = huffman.getHuffmanCode();

    // Step 6: Write frequencies and symbol count
    std::ofstream outFile(outputFile, std::ios::binary);
    if (!outFile) {
        std::cerr << "Cannot open output file: " << outputFile << std::endl;
        return;
    }

    // Write number of symbols (uint64_t)
    uint64_t numSymbols = mtfEncoded.size();
    outFile.write(reinterpret_cast<const char *>(&numSymbols), sizeof(numSymbols));

    // Write primary index for inverse BWT
    outFile.write(reinterpret_cast<const char *>(&primary_index), sizeof(primary_index));

    // Write frequencies
    for (int i = 0; i < 256; ++i) {
        uint64_t frequency = freq[i];
        outFile.write(reinterpret_cast<const char *>(&frequency), sizeof(frequency));
    }

    // Step 7: Write Huffman encoded symbols
    BitWriter writer(outFile);
    for (int symbol : mtfEncoded) {
        const auto &code = huffmanCode.at(symbol);
        for (bool bit : code) {
            writer.writeBit(bit);
        }
    }
    writer.flush(); // Ensure all bits are written

    outFile.close();
}

// Decompression Function
void decompress(const std::string &compressedFile, const std::string &outputFile) {
    std::ifstream inFile(compressedFile, std::ios::binary);
    if (!inFile) {
        std::cerr << "Cannot open compressed file: " << compressedFile << std::endl;
        return;
    }

    // Read number of symbols
    uint64_t numSymbols;
    inFile.read(reinterpret_cast<char *>(&numSymbols), sizeof(numSymbols));

    // Read primary index for inverse BWT
    int primary_index;
    inFile.read(reinterpret_cast<char *>(&primary_index), sizeof(primary_index));

    // Read frequencies
    std::map<int, uint64_t> freq;
    for (int i = 0; i < 256; ++i) {
        uint64_t frequency;
        inFile.read(reinterpret_cast<char *>(&frequency), sizeof(frequency));
        freq[i] = frequency;
    }

    // Build Huffman tree
    HuffmanDecoder huffmanDecoder(freq);

    // Use BitReader starting from the current position
    BitReader reader(inFile);

    // Decode Huffman data
    std::vector<int> mtfDecoded;
    try {
        mtfDecoded = huffmanDecoder.decode(reader, numSymbols);
    } catch (const std::runtime_error &e) {
        std::cerr << "Decompression error: " << e.what() << std::endl;
        return;
    }

    inFile.close();

    // MTF Decoding
    std::string bwtDecoded = move_to_front_decode(mtfDecoded);

    // Inverse BWT
    std::string original = inverse_burrows_wheeler_transform(bwtDecoded, primary_index);

    // Write the decompressed data to the output file
    std::ofstream outFile(outputFile, std::ios::binary);
    if (!outFile) {
        std::cerr << "Cannot open output file: " << outputFile << std::endl;
        return;
    }
    outFile.write(original.data(), original.size());
    outFile.close();
}

int main() {
    // Compress the input file
    compress("input.txt", "compressed.dat");

    // Decompress the compressed file
    decompress("compressed.dat", "decompressed.txt");

    std::cout << "Compression and decompression of 'input.txt' completed." << std::endl;

    // Compress the enwik8 file
    compress("enwik8.txt", "compressedenwik8.dat");

    // Decompress the compressed enwik8 file
    decompress("compressedenwik8.dat", "decompressedenwik8.txt");

    std::cout << "Compression and decompression of 'enwik8.txt' completed." << std::endl;
    return 0;
}
