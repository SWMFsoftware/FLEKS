#ifndef _BIT_H_
#define _BIT_H_

namespace bit {
//=== The bit that repsents the status of a cell/node. ===
constexpr static int iBny_ = 0;
constexpr static int iNew_ = 1;
constexpr static int iEdge_ = 2;
//=========================================================

inline bool test_bit(const int& i, int pos) { return i & (1 << pos); }

inline void turn_on_bit(int& i, int pos) { i |= (1 << pos); }

inline void turn_off_bit(int& i, int pos) { i &= ~(1 << pos); }

// inline void flip_bit(int& i, int pos) { i ^= (1 << pos); }

inline void set_boundary(int& i) { turn_on_bit(i, iBny_); }
inline void set_not_boundary(int& i) { turn_off_bit(i, iBny_); }
inline bool is_boundary(const int& i) { return test_bit(i, iBny_); }

inline void set_new(int& i) { turn_on_bit(i, iNew_); }
inline void set_not_new(int& i) { turn_off_bit(i, iNew_); }
inline bool is_new(const int& i) { return test_bit(i, iNew_); }

inline void set_edge(int& i) { turn_on_bit(i, iEdge_); }
inline void set_not_edge(int& i) { turn_off_bit(i, iEdge_); }
inline bool is_edge(const int& i) { return test_bit(i, iEdge_); }

} // namespace bit
#endif