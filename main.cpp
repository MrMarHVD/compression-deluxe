#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <queue>
#include <unordered_map>
#include <sstream>
#include <iterator>
#include <numeric>

using namespace std;

// Function to find a suitable EOF character
char findEOFChar(const string& data) {
    vector<bool> charExists(256, false);
    for (char c : data) {
        charExists[(unsigned char)c] = true;
    }
    for (int i = 1; i < 256; ++i) { // Start from 1 to avoid using '\0'
        if (!charExists[i]) {
            return (char)i;
        }
    }
    cerr << "No suitable EOF character found!" << endl;
    exit(1);
}

// Structure to store the BWT result
struct BWTResultat {
    string transformertData;
    int rotasjonsIndeks;
};

// Function to read the content of a file as binary data
string lesFil(const string& filnavn) {
    ifstream fil(filnavn, ios::binary);
    if (!fil) {
        cerr << "Could not open file: " << filnavn << endl;
        exit(1);
    }
    string innhold((istreambuf_iterator<char>(fil)), istreambuf_iterator<char>());
    return innhold;
}

// Function to write data to a file as binary data
void skrivFil(const string& filnavn, const string& data) {
    ofstream fil(filnavn, ios::binary);
    if (!fil) {
        cerr << "Could not write to file: " << filnavn << endl;
        exit(1);
    }
    fil.write(data.data(), data.size());
}

// Build suffix array
vector<int> buildSuffixArray(const string& s) {
    int n = s.size();
    vector<int> sa(n), ranks(n), tmp(n);

    // Initial ranking based on the first character
    for (int i = 0; i < n; ++i) {
        sa[i] = i;
        ranks[i] = (unsigned char)s[i];
    }

    for (int k = 1; k < n; k <<= 1) {
        auto cmp = [&](int i, int j) {
            if (ranks[i] != ranks[j])
                return ranks[i] < ranks[j];
            int ri = (i + k < n) ? ranks[i + k] : -1;
            int rj = (j + k < n) ? ranks[j + k] : -1;
            return ri < rj;
        };
        sort(sa.begin(), sa.end(), cmp);

        tmp[sa[0]] = 0;
        for (int i = 1; i < n; ++i) {
            tmp[sa[i]] = tmp[sa[i - 1]] + (cmp(sa[i - 1], sa[i]) ? 1 : 0);
        }
        ranks = tmp;
        if (ranks[sa[n - 1]] == n - 1)
            break; // All ranks are unique
    }
    return sa;
}

// Efficient Burrows-Wheeler Transform function
BWTResultat burrowsWheelerTransform(const string& input) {
    int n = input.size();
    vector<int> sa = buildSuffixArray(input);

    // Generate BWT from suffix array
    string bwt(n, '\0');
    int rotationIndex = -1;
    for (int i = 0; i < n; ++i) {
        if (sa[i] == 0)
            rotationIndex = i;
        bwt[i] = input[(sa[i] + n - 1) % n];
    }
    return {bwt, rotationIndex};
}

// Inverse Burrows-Wheeler Transform function
string inversBurrowsWheelerTransform(const string& bwt, int rotasjonsIndeks) {
    int n = bwt.size();
    vector<int> count(256, 0);

    // Count occurrence of each character in bwt
    for (int i = 0; i < n; ++i) {
        count[(unsigned char)bwt[i]]++;
    }

    // Compute cumulative counts
    vector<int> cumulativeCount(256, 0);
    for (int i = 1; i < 256; ++i) {
        cumulativeCount[i] = cumulativeCount[i - 1] + count[i - 1];
    }

    // Build rank array
    vector<int> rank(n, 0);
    vector<int> tally(256, 0);
    for (int i = 0; i < n; ++i) {
        rank[i] = tally[(unsigned char)bwt[i]];
        tally[(unsigned char)bwt[i]]++;
    }

    // Build LF mapping
    vector<int> LF(n, 0);
    for (int i = 0; i < n; ++i) {
        LF[i] = cumulativeCount[(unsigned char)bwt[i]] + rank[i];
    }

    // Reconstruct original string in reverse order
    string original(n, '\0');
    int idx = rotasjonsIndeks;
    for (int i = n - 1; i >= 0; --i) {
        original[i] = bwt[idx];
        idx = LF[idx];
    }

    return original;
}

// Move-to-Front Encoding function
vector<uint8_t> moveToFrontEncoding(const string& input) {
    vector<uint8_t> alfabet(256);
    iota(alfabet.begin(), alfabet.end(), 0);

    vector<uint8_t> output;
    output.reserve(input.size());

    for (unsigned char c : input) {
        uint8_t index = 0;
        while (alfabet[index] != c) {
            ++index;
        }
        output.push_back(index);

        // Move the found character to the front
        alfabet.erase(alfabet.begin() + index);
        alfabet.insert(alfabet.begin(), c);
    }

    return output;
}

// Move-to-Front Decoding function
string moveToFrontDecoding(const vector<uint8_t>& input) {
    vector<uint8_t> alfabet(256);
    iota(alfabet.begin(), alfabet.end(), 0);

    string output;
    output.reserve(input.size());

    for (uint8_t index : input) {
        unsigned char c = alfabet[index];
        output += c;

        // Move the found character to the front
        alfabet.erase(alfabet.begin() + index);
        alfabet.insert(alfabet.begin(), c);
    }

    return output;
}

// Node class for the Huffman tree
struct Node {
    uint8_t tegn;
    uint64_t frekvens;
    Node* venstre;
    Node* hoyre;

    Node(uint8_t t, uint64_t f) : tegn(t), frekvens(f), venstre(nullptr), hoyre(nullptr) {}
    Node(Node* v, Node* h) : tegn(0), frekvens(v->frekvens + h->frekvens), venstre(v), hoyre(h) {}
};

// Comparison function for the priority queue
struct NodeCompare {
    bool operator()(Node* a, Node* b) {
        return a->frekvens > b->frekvens;
    }
};

// Builds the Huffman tree from the frequency table
Node* byggHuffmanTre(const vector<uint64_t>& frekvenser) {
    priority_queue<Node*, vector<Node*>, NodeCompare> pq;

    for (uint16_t i = 0; i < 256; ++i) {
        if (frekvenser[i] > 0) {
            pq.push(new Node(static_cast<uint8_t>(i), frekvenser[i]));
        }
    }

    if (pq.empty()) {
        return nullptr;
    }

    while (pq.size() > 1) {
        Node* venstre = pq.top(); pq.pop();
        Node* hoyre = pq.top(); pq.pop();
        Node* forelder = new Node(venstre, hoyre);
        pq.push(forelder);
    }

    return pq.top();
}

// Builds the code table from the Huffman tree
void byggKoder(Node* rot, vector<string>& koder, string kode) {
    if (!rot->venstre && !rot->hoyre) {
        koder[rot->tegn] = kode;
        return;
    }
    byggKoder(rot->venstre, koder, kode + "0");
    byggKoder(rot->hoyre, koder, kode + "1");
}

// Frees the memory used by the Huffman tree
void slettTre(Node* rot) {
    if (!rot) return;
    slettTre(rot->venstre);
    slettTre(rot->hoyre);
    delete rot;
}

// Performs Huffman encoding
vector<uint8_t> huffmanKoding(const vector<uint8_t>& input, vector<string>& koder, string& frekvensData, uint32_t& bitLength) {
    // Calculate frequencies
    vector<uint64_t> frekvenser(256, 0);
    for (uint8_t c : input) {
        frekvenser[c]++;
    }

    // Build Huffman tree
    Node* rot = byggHuffmanTre(frekvenser);

    // Handle empty or single character case
    if (!rot) {
        return {};
    }

    // Build code table
    koder.resize(256);
    byggKoder(rot, koder, "");

    // Encode the data
    string kodetStreng;
    kodetStreng.reserve(input.size()); // To avoid frequent reallocations
    for (uint8_t c : input) {
        kodetStreng += koder[c];
    }

    // Store the total number of bits
    bitLength = kodetStreng.size();

    // Create frequency data for storage
    stringstream ss;
    for (uint16_t i = 0; i < 256; ++i) {
        if (frekvenser[i] > 0) {
            ss << i << ' ' << frekvenser[i] << ' ';
        }
    }
    frekvensData = ss.str();

    // Convert bit string to bytes
    vector<uint8_t> komprimertData;
    komprimertData.reserve((kodetStreng.size() + 7) / 8);
    uint8_t byte = 0;
    int bitsFyllt = 0;
    for (char bit : kodetStreng) {
        byte = (byte << 1) | (bit - '0');
        bitsFyllt++;
        if (bitsFyllt == 8) {
            komprimertData.push_back(byte);
            bitsFyllt = 0;
            byte = 0;
        }
    }
    if (bitsFyllt > 0) {
        byte <<= (8 - bitsFyllt);
        komprimertData.push_back(byte);
    }

    slettTre(rot);
    return komprimertData;
}

// Performs Huffman decoding
vector<uint8_t> huffmanDekoding(const vector<uint8_t>& kodetData, vector<string>& koder, const string& frekvensData, uint32_t bitLength) {
    // Rebuild the frequency table
    vector<uint64_t> frekvenser(256, 0);
    stringstream ss(frekvensData);
    uint16_t tegn;
    uint64_t frekvens;
    while (ss >> tegn >> frekvens) {
        frekvenser[tegn] = frekvens;
    }

    // Build Huffman tree
    Node* rot = byggHuffmanTre(frekvenser);

    // Handle empty or single character case
    if (!rot) {
        return {};
    }

    // Build decode table
    unordered_map<string, uint8_t> dekodeTabell;
    koder.resize(256);
    byggKoder(rot, koder, "");
    for (uint16_t i = 0; i < 256; ++i) {
        if (!koder[i].empty()) {
            dekodeTabell[koder[i]] = static_cast<uint8_t>(i);
        }
    }

    // Reconstruct the bitstring up to the valid bit length
    string bitstring;
    bitstring.reserve(bitLength);
    int bitsRead = 0;
    for (uint8_t byte : kodetData) {
        for (int i = 7; i >= 0 && bitsRead < bitLength; --i) {
            bitstring += ((byte >> i) & 1) ? '1' : '0';
            bitsRead++;
        }
    }

    vector<uint8_t> output;
    output.reserve(bitLength / 8); // Approximate size
    string kode;
    for (char bit : bitstring) {
        kode += bit;
        if (dekodeTabell.count(kode)) {
            output.push_back(dekodeTabell[kode]);
            kode.clear();
        }
    }

    slettTre(rot);
    return output;
}

// Compresses the input file and saves it as the output file
void komprimer(const string& inputFilnavn, const string& outputFilnavn) {
    // Read data from file
    string data = lesFil(inputFilnavn);

    // Find EOF character
    char EOF_char = findEOFChar(data);

    // Append EOF character
    data += EOF_char;

    // Perform BWT
    BWTResultat bwtResultat = burrowsWheelerTransform(data);

    // Perform Move-to-Front Encoding
    vector<uint8_t> mtfData = moveToFrontEncoding(bwtResultat.transformertData);

    // Perform Huffman encoding
    vector<string> koder;
    string frekvensData;
    uint32_t bitLength; // Variable to store the bit length
    vector<uint8_t> komprimertData = huffmanKoding(mtfData, koder, frekvensData, bitLength);

    // Save rotation index, frequency data, bit length, EOF character, and compressed data to file
    ofstream fil(outputFilnavn, ios::binary);
    if (!fil) {
        cerr << "Could not write to file: " << outputFilnavn << endl;
        exit(1);
    }
    // Write rotation index
    fil.write(reinterpret_cast<const char*>(&bwtResultat.rotasjonsIndeks), sizeof(bwtResultat.rotasjonsIndeks));

    // Write length of frequency data
    uint32_t frekvensDataLengde = frekvensData.size();
    fil.write(reinterpret_cast<const char*>(&frekvensDataLengde), sizeof(frekvensDataLengde));

    // Write frequency data
    fil.write(frekvensData.data(), frekvensData.size());

    // Write bit length
    fil.write(reinterpret_cast<const char*>(&bitLength), sizeof(bitLength));

    // Write EOF character
    fil.write(&EOF_char, sizeof(EOF_char));

    // Write compressed data
    fil.write(reinterpret_cast<const char*>(komprimertData.data()), komprimertData.size());

    fil.close();
}

// Decompresses the input file and recreates the output file
void dekomprimer(const string& inputFilnavn, const string& outputFilnavn) {
    // Read compressed data from file
    ifstream fil(inputFilnavn, ios::binary);
    if (!fil) {
        cerr << "Could not open file: " << inputFilnavn << endl;
        exit(1);
    }

    // Read rotation index
    int rotasjonsIndeks;
    fil.read(reinterpret_cast<char*>(&rotasjonsIndeks), sizeof(rotasjonsIndeks));

    // Read length of frequency data
    uint32_t frekvensDataLengde;
    fil.read(reinterpret_cast<char*>(&frekvensDataLengde), sizeof(frekvensDataLengde));

    // Read frequency data
    string frekvensData(frekvensDataLengde, '\0');
    fil.read(&frekvensData[0], frekvensDataLengde);

    // Read the bit length
    uint32_t bitLength;
    fil.read(reinterpret_cast<char*>(&bitLength), sizeof(bitLength));

    // Read EOF character
    char EOF_char;
    fil.read(&EOF_char, sizeof(EOF_char));

    // Read the rest as compressed data
    vector<uint8_t> komprimertData((istreambuf_iterator<char>(fil)), istreambuf_iterator<char>());
    fil.close();

    // Perform Huffman decoding
    vector<string> koder;
    vector<uint8_t> mtfData = huffmanDekoding(komprimertData, koder, frekvensData, bitLength);

    // Perform Move-to-Front Decoding
    string bwtData = moveToFrontDecoding(mtfData);

    // Perform inverse BWT
    string originalData = inversBurrowsWheelerTransform(bwtData, rotasjonsIndeks);

    // Remove EOF character from the end
    if (!originalData.empty() && originalData.back() == EOF_char) {
        originalData.pop_back();
    }

    // Write decompressed data to file
    skrivFil(outputFilnavn, originalData);
}

int main() {
    string inputFil = "input.txt";
    string komprimertFil = "komprimert.dat";
    string dekomprimertFil = "dekomprimert.txt";

    komprimer(inputFil, komprimertFil);
    dekomprimer(komprimertFil, dekomprimertFil);


    string inputFil2 = "diverse.lyx";
    string komprimertFil2 = "komprimert2.dat";
    string dekomprimertFil2 = "dekomprimert2.lyx";

    komprimer(inputFil2, komprimertFil2);
    dekomprimer(komprimertFil2, dekomprimertFil2);


    return 0;
}
