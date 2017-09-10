#include <cstddef>

#include <gmpxx.h>
#include <algorithm>
#include <chrono>
#include <nfl.hpp>
#include <time.h>
#include "utils.h"

using namespace std;

/// include the FV homomorphic encryption library
namespace FV {
  namespace params {
    using poly_t = nfl::poly_from_modulus<uint64_t, 1 << 12, 248>;
    
    template <typename T>
    struct plaintextModulus;
    
    template <>
    struct plaintextModulus<mpz_class> 
    {
      static mpz_class value() 
      { 
        return mpz_class(2);
      }
    };
    
    template <>
    struct plaintextModulus<unsigned long> 
    {
      static unsigned long value() 
      {
        return 2;
      }
    };
    
    using gauss_struct = nfl::gaussian<uint16_t, uint64_t, 2>;
    using gauss_t = nfl::FastGaussianNoise<uint16_t, uint64_t, 2>;
    
    gauss_t fg_prng_sk(8.0, 128, 1 << 14);
    gauss_t fg_prng_evk(8.0, 128, 1 << 14);
    gauss_t fg_prng_pk(8.0, 128, 1 << 14);
    gauss_t fg_prng_enc(8.0, 128, 1 << 14);
  }
}  // namespace FV::params
#include <FV.hpp>

using namespace FV;

const int OSIZE = 10;
const int N1 = 42;
const int N2 = 128;
const int N3 = 9;
const int N4 = 8;
const int N3_sum = 45;
const int SIZE = N1 + N2 + (N3_sum * N4);

void permutation(std::array<int, SIZE> &values){  
  int nb, tmp;
  for ( int i = 0; i < SIZE; i++ ) values[i] = i;
  
  for (int i = SIZE - 1; i > 0; i--){
    nb = rand() & i;
    tmp = values[i];
    values[i] = values[nb];
    values[nb] = tmp;
  }
}

template<typename S>
S linear(std::array<S, SIZE> const &state){  
  S result = state[0];

  for (int i = 1; i < N1; i++) result +=state[i];

  return result;
}

template<typename S>
S quadratic(std::array<S, SIZE> const &state){
  S result = state[N1] * state[N1 + 1];
  for (int i = 2; i < N2; i += 2) result += state[N1 + i] * state[N1 + i + 1];

  return result;
}

template<typename S>
S triangular(std::array<S, SIZE> const &state){
  int it = N1 + N2;
  S result = state[it];

  for(int i = 0; i < N4; i++){
    int multiplications = 0;

    for(it = N1 + N2 + (i * N3_sum); it < N1 + N2 + ((i+1) * N3_sum); ){
      multiplications+= 1;
      S tmp = state[it];
      //cout << it << " ";
      it++;
      for(int j = 0; j < multiplications - 1; j++){
        tmp *= state[it];
        //cout << it << " ";
        it++;
      }
      //cout << endl;
      result += tmp;

    }
  }

  return result;
}

template<typename K, typename S, typename O>
void Flip(std::array<K, SIZE> const &key, std::array<S, SIZE> &state, std::array<O, OSIZE> &output)
{

  std::array<int, SIZE> pm;

  for (int i = 0; i < SIZE; i++) state[i] = key[i];

  std::array<S, SIZE> alter_state;
  
  //Iterate until you have enough output
  for (int i = 0; i < OSIZE; i++){

    //Ask for new permutation
    permutation(pm);

    /*for (auto const &v : pm) {
      std::cout << v << " ";
    }
    cout << endl;*/

    //Copy initial state to altered state using generated permutation
    for (int j = 0; j < SIZE; j++) alter_state[pm[j]] = state[j];

    output[i] = linear(alter_state) + quadratic(alter_state) + triangular(alter_state);
  }

}

int main()
{                     
  FV::sk_t sk;
  FV::evk_t evk(sk, 32);
  FV::pk_t pk(sk, evk);

  // Timer
  auto start = std::chrono::steady_clock::now();
  auto end = std::chrono::steady_clock::now();
  
  srand(time(NULL));
  std::array<unsigned long, SIZE> key; 
  std::generate(key.begin(), key.end(), [] { return rand() & 1; });
  
  std::array<FV::message_t<unsigned long>, SIZE> state;
  std::array<FV::message_t<unsigned long>, OSIZE> output;
  
  srand(0);
  start = std::chrono::steady_clock::now();
  Flip(key, state, output);
  end = std::chrono::steady_clock::now();

  std::cout << "\tBinary Flip: \t\t\t"
            << get_time_us(start, end, 2) << " us"
            << std::endl;


  for (auto const &v : output) {
    std::cout << v;
  }
  cout << endl;
 
  std::array<FV::ciphertext_t, SIZE> h_key;
  for (size_t i = 0; i < SIZE; i++) {
  if (key[i] == 1)
    encrypt_integer(h_key[i], pk, key[i]);
  else
    h_key[i] = 0;
  }

  std::array<FV::ciphertext_t, SIZE> h_state;
  std::generate(h_state.begin(), h_state.end(), [&pk] {
    return FV::ciphertext_t(pk, 0);
  });  // set the public key in all the state ciphertexts

  std::array<FV::ciphertext_t, OSIZE> h_output;

  srand(0);
  start = std::chrono::steady_clock::now();
  Flip(h_key, h_state, h_output);
  end = std::chrono::steady_clock::now();

  std::cout << "\tHomomorphic Flip: \t\t"
            << get_time_us(start, end, 2) / 1000000 << " s"
            << std::endl;

  std::array<FV::mess_t, OSIZE> d_output;
  for (int i = 0; i < OSIZE; i++){
    //h_output[i]+=1;
    decrypt(d_output[i], sk, pk, h_output[i]);
    std::cout << d_output[i];
  }
  cout << endl;

  return 0;     
}
