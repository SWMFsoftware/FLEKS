#ifndef _BIT_H_
#define _BIT_H_

namespace bit {
//=== The bit that repsents the status of a cell/node. ===

// Example: iBny_ = 0. It means the 0th bit is used to represent if a cell/node
// is a boundary cell/node. If the 0th bit is 1, then it is a boundary
// cell/node.

// Boundary cell/node. Note: the nodes just at the boundary are not 'boundary'
constexpr static int iBny_ = 0;

// New active cell/node
constexpr static int iNew_ = 1;

// The nodes at the boundary or the physical cells just next to the boundary.
constexpr static int iEdge_ = 2;

// Node only. A node may be shared by a few blocks/boxes. Sometimes (such as the
// E field solver) only one of the boexes needs to do some operations for this
// node. This box is the 'owner' of this node.
constexpr static int iOwner_ = 3;

// For cell only. If a cell is refined or not.
constexpr static int iRefined_ = 4;
//=========================================================

inline bool test_bit(const int& i, int pos) { return i & (1 << pos); }
inline void turn_on_bit(int& i, int pos) { i |= (1 << pos); }
inline void turn_off_bit(int& i, int pos) { i &= ~(1 << pos); }

// inline void flip_bit(int& i, int pos) { i ^= (1 << pos); }

//======= Boundary =======
inline void set_boundary(int& i) { turn_on_bit(i, iBny_); }
inline void set_not_boundary(int& i) { turn_off_bit(i, iBny_); }
inline bool is_boundary(const int& i) { return test_bit(i, iBny_); }

//======= New active cell/node =======
inline void set_new(int& i) { turn_on_bit(i, iNew_); }
inline void set_not_new(int& i) { turn_off_bit(i, iNew_); }
inline bool is_new(const int& i) { return test_bit(i, iNew_); }

//======= Edge cell/node =======
inline void set_edge(int& i) { turn_on_bit(i, iEdge_); }
inline void set_not_edge(int& i) { turn_off_bit(i, iEdge_); }
inline bool is_edge(const int& i) { return test_bit(i, iEdge_); }

//======= Owner node (node only) =======
inline void set_owner(int& i) { turn_on_bit(i, iOwner_); }
inline void set_not_owner(int& i) { turn_off_bit(i, iOwner_); }
inline bool is_owner(const int& i) { return test_bit(i, iOwner_); }

//======= Refined cell =======
inline void set_refined(int& i) { turn_on_bit(i, iRefined_); }
inline void set_not_refined(int& i) { turn_off_bit(i, iRefined_); }
inline bool is_refined(const int& i) { return test_bit(i, iRefined_); }

} // namespace bit
#endif