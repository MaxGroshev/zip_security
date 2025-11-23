#include <iostream>
#include <vector>
#include <array>
#include <cstdint>
#include <cstring>
#include <iomanip>

class ChaCha20 {
public:
    // ChaCha20 constants
    static const int KEY_SIZE = 32;
    static const int NONCE_SIZE = 12; // 96-bit nonce (RFC 8439)
    static const int BLOCK_SIZE = 64;

    /**
     * @brief Initialize ChaCha20 state
     * @param key 32-byte key
     * @param nonce 12-byte nonce
     * @param counter Initial counter (usually 0 or 1)
     */
    ChaCha20(const std::array<uint32_t, KEY_SIZE>& key, 
             const std::array<uint32_t, NONCE_SIZE>& nonce, 
             uint32_t counter = 1) {
        stateSetup(key, nonce, counter);
    }

    /**
     * @brief Encrypt or Decrypt data (XOR operation is symmetric)
     * @param input Input data buffer
     * @param output Output data buffer (resized automatically)
     */
    void crypt(const std::vector<uint32_t>& input, std::vector<uint32_t>& output) {
        output.resize(input.size());
        std::vector<uint32_t> keyStreamBlock(BLOCK_SIZE);
        
        // Process full blocks
        size_t numBlocks = input.size() / BLOCK_SIZE;
        size_t remainder = input.size() % BLOCK_SIZE;

        for (size_t i = 0; i < numBlocks; ++i) {
            generateBlock(keyStreamBlock);
            for (size_t j = 0; j < BLOCK_SIZE; ++j) {
                output[i * BLOCK_SIZE + j] = input[i * BLOCK_SIZE + j] ^ keyStreamBlock[j];
            }
        }

        // Process remaining bytes
        if (remainder > 0) {
            generateBlock(keyStreamBlock);
            for (size_t j = 0; j < remainder; ++j) {
                output[numBlocks * BLOCK_SIZE + j] = input[numBlocks * BLOCK_SIZE + j] ^ keyStreamBlock[j];
            }
        }
    }

private:
    uint32_t state[16];

    // Bit rotation helper
    static inline uint32_t rotl(uint32_t x, int n) {
        return (x << n) | (x >> (32 - n));
    }

    // Pack 4 bytes into a 32-bit integer (Little Endian)
    static inline uint32_t pack4(const uint32_t* a) {
        return uint32_t(a[0]) | (uint32_t(a[1]) << 8) | 
               (uint32_t(a[2]) << 16) | (uint32_t(a[3]) << 24);
    }

    // Unpack 32-bit integer into 4 bytes (Little Endian)
    static inline void unpack4(uint32_t src, std::vector<uint32_t>& dst, size_t offset) {
        dst[offset + 0] = (src >> 0) & 0xff;
        dst[offset + 1] = (src >> 8) & 0xff;
        dst[offset + 2] = (src >> 16) & 0xff;
        dst[offset + 3] = (src >> 24) & 0xff;
    }

    void stateSetup(const std::array<uint32_t, KEY_SIZE>& key, 
                    const std::array<uint32_t, NONCE_SIZE>& nonce, 
                    uint32_t counter) {
        // 1. Constants "expand 32-byte k"
        state[0] = 0x61707865;
        state[1] = 0x3320646e;
        state[2] = 0x79622d32;
        state[3] = 0x6b206574;

        // 2. Key (32 bytes -> 8 words)
        for (int i = 0; i < 8; ++i) {
            state[4 + i] = pack4(&key[i * 4]);
        }

        // 3. Counter (32-bit)
        state[12] = counter;

        // 4. Nonce (12 bytes -> 3 words)
        for (int i = 0; i < 3; ++i) {
            state[13 + i] = pack4(&nonce[i * 4]);
        }
    }

    // The Quarter Round function
    static void quarterRound(uint32_t& a, uint32_t& b, uint32_t& c, uint32_t& d) {
        a += b; d ^= a; d = rotl(d, 16);
        c += d; b ^= c; b = rotl(b, 12);
        a += b; d ^= a; d = rotl(d, 8);
        c += d; b ^= c; b = rotl(b, 7);
    }

    void generateBlock(std::vector<uint32_t>& outputBlock) {
        uint32_t workingState[16];
        std::memcpy(workingState, state, sizeof(state));

        // 10 double-rounds (20 rounds total)
        for (int i = 0; i < 10; ++i) {
            // Column rounds
            quarterRound(workingState[0], workingState[4], workingState[8],  workingState[12]);
            quarterRound(workingState[1], workingState[5], workingState[9],  workingState[13]);
            quarterRound(workingState[2], workingState[6], workingState[10], workingState[14]);
            quarterRound(workingState[3], workingState[7], workingState[11], workingState[15]);

            // Diagonal rounds
            quarterRound(workingState[0], workingState[5], workingState[10], workingState[15]);
            quarterRound(workingState[1], workingState[6], workingState[11], workingState[12]);
            quarterRound(workingState[2], workingState[7], workingState[8],  workingState[13]);
            quarterRound(workingState[3], workingState[4], workingState[9],  workingState[14]);
        }

        // Add original state to working state
        for (int i = 0; i < 16; ++i) {
            workingState[i] += state[i];
        }

        // Serialize (Little Endian) to output stream
        for (int i = 0; i < 16; ++i) {
            unpack4(workingState[i], outputBlock, i * 4);
        }

        // Increment block counter
        state[12]++;
        // Note: If state[12] overflows, it wraps to 0.
        // RFC 8439 doesn't specify using the nonce as the upper 32 bits of a 64-bit counter,
        // but some implementations stop if the counter limits are reached.
    }
};

// --- Helper for printing ---
void printHex(const std::vector<uint32_t>& data) {
    for(auto byte : data) {
        std::cout << std::hex << std::setw(2) << std::setfill('0') << (int)byte;
    }
    std::cout << std::dec << std::endl;
}

int main() {
    // Example Test Vector from RFC 8439 (Sunscreen)
    
    // Key: 00:01:02:03:04:05:06:07:08:09:0a:0b:0c:0d:0e:0f...1f
    std::array<uint32_t, 32> key;
    for(int i=0; i<32; i++) key[i] = i+1;

    // Nonce: 00:00:00:00:00:00:00:4a:00:00:00:00
    std::array<uint32_t, 12> nonce = {0,0,0,0,0,0,0,0x4a,0,0,0,0};
    
    // Initial counter
    uint32_t counter = 1;

    // Message: "Ladies and gentlemen of the class of '99: If I could offer you only one tip for the future, sunscreen would be it."
    std::string msgStr = "Ladies and gentlemen of the class of '99: If I could offer you only one tip for the future, sunscreen would be it.";
    std::vector<uint32_t> plaintext(msgStr.begin(), msgStr.end());
    
    std::cout << "Original Text: " << msgStr << "\n\n";

    // 1. Encryption
    ChaCha20 chacha(key, nonce, counter);
    std::vector<uint32_t> ciphertext;
    chacha.crypt(plaintext, ciphertext);

    std::cout << "Ciphertext (Hex): ";
    printHex(ciphertext);

    // 2. Decryption (Create new instance with same Key/Nonce/Counter)
    ChaCha20 chachaDecryptor(key, nonce, counter);
    std::vector<uint32_t> decryptedText;
    chachaDecryptor.crypt(ciphertext, decryptedText);

    std::string recoveredStr(decryptedText.begin(), decryptedText.end());
    std::cout << "\nDecrypted Text: " << recoveredStr << "\n";
    
    return 0;
}