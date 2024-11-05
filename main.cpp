#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <queue>
#include <bitset>
#include <cstdint>
#include <string>


// Node of the Huffman tree
// Leafs represent characters
// Internal nodes represent frequencies
// Nodes are organized by frequency when doing Huffman compression
class TreeNode {
    public:
        char byte;
        int frequency;
        TreeNode *left; // Pointer to the node's left child
        TreeNode *right; // Pointer to the node's right child
        TreeNode(char b, int freq) : byte(b), frequency(freq), left(nullptr), right(nullptr) {}
};
// Used to compare nodes in the priority queue when generating a Huffman tree
// Priority should be given to low frequency nodes
class Compare {
    public:
    bool operator()(TreeNode* left, TreeNode* right) {
        return left-> frequency > right-> frequency;
    }
};

// Compression Class based on Huffman Encoding
// Huffman Coding is the second step of the compression process
// The Compression class builds a Huffman tree based on character frequency

class Compression {

    public:
        TreeNode* root; // Root of the Huffman tree

    Compression() : root(nullptr) {}
    // Serializes Huffman tree
    // '1' represents a non-null node
    // '0' represents a null node
    //TODO Replace characters with something more efficient?
    void serializeTree(TreeNode* node, std::ostream& output) {
        if (!node) {
            output.put('0'); // Empty child
            return;
        }
        output.put('1'); // Non-null child
        output.put(node->byte); // Stores byte value at leaf node
        serializeTree(node->left, output); // Recursive serialization of left subtree
        serializeTree(node->right, output);  // Recursive serialization of right subtree
    }

    // Deserializes Huffman tree for decompression
    // Reads binary format of tree and reconstructs it
    //TODO Replace characters with something more efficient?
    TreeNode* deserializeTree(std::ifstream& input) {
        char marker;
        input.get(marker);
        if (marker == '0') return nullptr; // Returns nullptr for null node
        char byte;
        input.get(byte); // Reads the byte value of the node
        TreeNode* node = new TreeNode(byte, 0); // Disregard frequency for decompression
        node->left = deserializeTree(input); // Reconstructs left subtree
        node->right = deserializeTree(input); // Reconstructs right subtree
        return node;
    }
    // Compresses input file with Huffman coding.
    // Compressed data is written to "compressed.txt"
    // TODO Make different compressed files for different input files
    void compressFile(const std::string &inputFile, const std::string &outputFile) {
        std::ifstream input(inputFile, std::ios::binary );
        std::map<char, int> frequencyTable;
        char byte;

        // Step 1: Build frequency table for bytes in file
        while (input.get(byte)) {
            frequencyTable[byte]++;
        }
        input.close();

        // Step 2: Construct Huffman tree from byte frequencies
        std::priority_queue<TreeNode*, std::vector<TreeNode*>, Compare> pq;
        for (auto& pair : frequencyTable) {
            pq.push(new TreeNode(pair.first, pair.second));
        }

        // Step 3: Combine lowest frequency nodes to from tree
        while (pq.size()>1) {
            TreeNode* left = pq.top();
            pq.pop();
            TreeNode* right = pq.top();
            pq.pop();
            TreeNode* combined = new TreeNode('\0', left-> frequency + right->frequency);
            combined->left = left;
            combined->right = right;
            pq.push(combined);
        }
        if (!pq.empty()) {
            root=pq.top(); // Root of Huffman tree
        }

        // Step 4: Generates Huffman codes for each byte
        std::map<char, std::string> huffmanCodes;
        generateCodes(root, "", huffmanCodes);

        // Step 5: Re-open input file to apply Huffman encoding
        input.open(inputFile, std::ios::binary);
        std::ofstream output(outputFile, std::ios::binary);

        // Step 5: Serialize Huffman tree structure
        serializeTree(root, output);

        // Step 6: Encode the content using Huffman codes
        // Stores content in bitstring
        std::string bitString;
        while (input.get(byte)) {
            bitString += huffmanCodes[byte];
        }
        input.close();

        // Step 7: Write the length of useful bits into the encodes bit string
        uint32_t bitLength = bitString.size();
        output.write(reinterpret_cast<char*>(&bitLength), sizeof(bitLength));

        // Step 8: Write encoded bitstring to file. Convert bits to bytes
        for (size_t i = 0; i < bitString.size(); i += 8) {
            std::bitset<8> bits(bitString.substr(i, 8).append(8 - (bitString.size() % 8), '0'));
            output.put(static_cast<unsigned char>(bits.to_ulong()));
        }
        output.close();
    }

    // Decompress Humman-compressd file
    // Store output in "Tempcompressed.txt"
    //TODO Create different file for different input filename
    void decompressFile(const std::string &inputFile) {
        std::ifstream input(inputFile, std::ios::binary);
        std::ofstream output("Tempcompressed.txt", std::ios::binary);

        // Step 1: Deserialize Huffman tree from compressed file
        root = deserializeTree(input);
        if (!root) {
            std::cerr << "Decompression error: Huffman tree not built" << std::endl;
            return;
        }

        // Step 2: Read the length of useful bits for accurate decompression
        uint32_t bitLength;
        input.read(reinterpret_cast<char*>(&bitLength), sizeof(bitLength));

        // Step 3: Rad compressed bit data from file and store it as bitstring
        std::string bitString;
        char byte;
        while (input.get(byte)) {
            bitString += std::bitset<8>(byte).to_string();
        }
        bitString = bitString.substr(0, bitLength);

        // Step 4: Decompress data using Huffman tree
        // Traverses tree bit by bit
        TreeNode* node = root;
        for (char bit : bitString) {
            node = (bit == '1') ? node->right : node->left;
            if (node && !node->left && !node->right) {
                output.put(node->byte); // Outputs byte when reaching a leaf node
                 node = root;
            }
        }
        input.close();
        output.close();
    }

    private:
    // Recursive helper function which generates Huffman codes for characters
    void generateCodes(TreeNode* node, const std::string& prefix, std::map<char, std::string>& codes) {
        if (!node) return;
        if (!node->left && !node->right) codes[node->byte] = prefix; //Stores code

        generateCodes(node->left, prefix + "0", codes);
        generateCodes(node->right, prefix + "1", codes);
    }
};

// Lempel-Ziv Compression Class
// Class handles compression and decompression using LZ algorithm
// Finds repeated byte sequences in a sliding window and replaces them with references
class LempelZivCompression {
    public:
        void initFiles(const std::string &inputFileName, const std::string & outputFileName) {
            inputFile.open(inputFileName, std::ios::binary);
            outputFile.open(outputFileName, std::ios::binary);
            if (!inputFile.is_open() || !outputFile.is_open()) {
                std::cerr << "Failed to open file" << std::endl;
                exit(1);
            }
        }

    // Compress data by identifying repeated sequences within a sliding window
    //TODO find best bufferSize - Possibly allow user choice for different file sizes
    void LempelZiv(size_t bufferSize) {
            std::vector<char> data;
            char byte;
            while (inputFile.get(byte)) {
                data.push_back(byte);
            }
            size_t dataSize = data.size();


            // Process each byte in input data
            for (size_t i = 0; i < dataSize; ++i) {
                size_t bestOffset = 0;
                size_t bestLength = 0;

                // Search for best match in buffer
                for (size_t j = 1; j <= bufferSize && j <= i; ++j) {
                    size_t length = 0;
                    while (length < bufferSize && i + length < dataSize && data[i - j + length] == data[i + length]) {
                        ++length;
                    }
                    if (length > bestLength && (i - j + length) <= i) {
                        bestOffset = j;
                        bestLength = length;
                    }
                }
                if (bestLength >= 4 && i >= bestOffset) {
                    // Writes reference to the match
                    uint16_t bestOffset16 = static_cast<uint16_t>(bestOffset);
                    uint16_t bestLength16 = static_cast<uint16_t>(bestLength);
                    outputFile.put(1); // Flag indicating match
                    outputFile.write(reinterpret_cast<const char*>(&bestOffset16), sizeof(bestOffset16));
                    outputFile.write(reinterpret_cast<const char*>(&bestLength16), sizeof(bestLength16));
                    i += bestLength - 1; // Skips past matched sequence


                } else {
                    // Write literal byte
                    outputFile.put(0); // Flag indicating literal byte
                    outputFile.put(data[i]);
                }

            }

        }
    void closeFiles() {
            inputFile.close();
            outputFile.close();
        }

    private:
        std::ifstream inputFile;
        std::ofstream outputFile;
};

// Lempel-Ziv Decompression
// Decompresses data using Lempel-Ziv algorithm
// Reconstructs sequences based on matched references and literal bytes
class LempelZivDecompression {
    public:
        void initFiles(const std::string &inputFileName, const std::string & outputFileName) {
            inputFile.open(inputFileName, std::ios::binary);
            outputFile.open(outputFileName, std::ios::binary);
            if (!inputFile.is_open()|| !outputFile.is_open()) {
                std::cerr << "Failed to open file" << std::endl;
                exit(1);
            }
        }

    void Decompression() {
            std::vector<char> data;
            uint16_t offset;
            uint16_t length;
            char literal;
            while (inputFile.peek() != EOF) {
                char flag;
                inputFile.get(flag);

                if (flag == 0) {
                    // Read and store literal byte
                    inputFile.get(literal);
                    data.push_back(literal);
                } else if (flag == 1) {
                    // Read offset and length
                    // Copy matched sequence
                    inputFile.read(reinterpret_cast<char*>(&offset), sizeof(offset));
                    inputFile.read(reinterpret_cast<char*>(&length), sizeof(length));
                    size_t start = data.size();
                    for (size_t i = 0; i < length; ++i) {
                        data.push_back(data[start - offset + i]);
                    }
                }
            }
            // Output for decompressed data
            for (char b : data) {
                outputFile.put(b);
            }
        }
    void closeFiles() {
            inputFile.close();
            outputFile.close();
        }
    private:
        std::ifstream inputFile;
        std::ofstream outputFile;
};

// Compression functions which encapsulate the main compressin steps
void LZCompressionStep(const std::string& fileName, size_t bufferSize) {
    LempelZivCompression lz;
    lz.initFiles(fileName, "Tempcompressed.txt");
    lz.LempelZiv(bufferSize);
    lz.closeFiles();
}
void LZDecompressionStep(const std::string& fileName) {
    LempelZivDecompression lz;
    lz.initFiles("Tempcompressed.txt", fileName);
    lz.Decompression();
    lz.closeFiles();
}
void HuffCompressionStep(const std::string& compressedFileName) {
    Compression huff;
    huff.compressFile("Tempcompressed.txt", compressedFileName);
}
void HuffDeCompressionStep(const std::string& compressedFileName) {
    Compression huff;
    huff.decompressFile(compressedFileName);
}
//TODO Better user interface
//TODO Split functionalit into two programs?
int main()
{
    std::string choice;
    std::cout << "The compression algorithm uses different buffer sizes for LZ for different files" << std::endl;
    std::cout << "This ensures better compression for diverse.lxy, but allows enwik8 to be compressed" << std::endl;
    std::cout << "Please choose an option:" << std::endl;
    std::cout << "1 - Compress"<< std::endl;
    std::cout << "2 - Decompress"<< std::endl;
    std::cout << "Hit another button to exit"<< std::endl;
    std::cin >> choice;
    if (choice == "1") {
        std::string choice2;
        std::cout << "Which file to compress?" << std::endl;
        std::cout << "1 - diverse.lyx"<< std::endl;
        std::cout << "2 - enwik8.txt"<< std::endl;
        std::cout << "3 - Twenty_thousand_leagues_under_the_sea.txt"<< std::endl;
        std::cout << "4 - opg6-kompr.lyx"<< std::endl;
        std::cout << "5 - diverse.txt"<< std::endl;



        std::cin >> choice2;
        if (choice2 == "1") {
            LZCompressionStep("diverse.lyx", 32768);
            HuffCompressionStep("diverse_compressed");
            std::cout << "File compressed to 'diverse_compressed'" << std::endl;
        } else if (choice2 == "2") {
            LZCompressionStep("enwik8.txt", 4096);
            HuffCompressionStep("enwik8compressed");
            std::cout << "File compressed to 'enwik8compressed'" << std::endl;
        } else if (choice2 == "3") {
            LZCompressionStep("Twenty_thousand_leagues_under_the_sea.txt", 32768);
            HuffCompressionStep("twenty_thousand_compressed");
            std::cout << "File compressed to 'twenty_thousand_compressed'" << std::endl;
        } else if (choice2 == "4") {
            LZCompressionStep("opg6-kompr.lyx", 4096);
            HuffCompressionStep("opg6-kompr_compressed");
            std::cout << "File compressed to 'opg6-kompr_compressed'" << std::endl;
        } else if (choice2 == "5") {
            LZCompressionStep("diverse.txt", 32768);
            HuffCompressionStep("diverse_txt_compressed");
            std::cout << "File compressed to 'diverse_txt_compressed'" << std::endl;
        }
    } else if (choice == "2") {
        std::string choice3;
        std::cout << "Which file to decompress?" << std::endl;
        std::cout << "1 - diverse.lyx"<< std::endl;
        std::cout << "2 - enwik8.txt"<< std::endl;
        std::cout << "3 - Twenty_thousand_leagues_under_the_sea.txt"<< std::endl;
        std::cout << "4 - opg6-kompr.lyx"<< std::endl;
        std::cout << "5 - diverse.txt"<< std::endl;
        std::cin >> choice3;
        if (choice3 == "1") {
            HuffDeCompressionStep("diverse_compressed");
            LZDecompressionStep("diverse_decompressed.lyx"
                            );
            std::cout << "File decompressed" << std::endl;
        } else if (choice3 == "2") {
            HuffDeCompressionStep("enwik8compressed");
            LZDecompressionStep("enwik8_decompressed.txt"
                            );
            std::cout << "File decompressed" << std::endl;
        } else if (choice3 == "3") {
            HuffDeCompressionStep("twenty_thousand_compressed");
            LZDecompressionStep("Twenty_thousand_leagues_under_the_sea_decompressed.txt"
                            );
            std::cout << "File decompressed" << std::endl;

        } else if (choice3 == "4") {
            HuffDeCompressionStep("opg6-kompr_compressed");
            LZDecompressionStep("opg6-kompr_sea_decompressed.txt"
                            );
            std::cout << "File decompressed" << std::endl;
        } else if (choice3 == "5") {
            HuffDeCompressionStep("diverse_txt_compressed");
            LZDecompressionStep("diverse_decompressed.txt"
                            );
            std::cout << "File decompressed" << std::endl;
        }



    } else {
        std::cout << "Exiting program" << std::endl;
    }
    return 0;
}
