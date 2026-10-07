# Useful sites: https://isocpp.org/wiki/faq

# Standards
1. Avoid using raw pointer and 'new' to allocate memory unless there is a good reason. 
   In other words, shared_ptr or unique_ptr should be always used to manage ownership. 
   If there is no notion of ownership, raw pointer should be used.

2. Naming standards:
 - File name: FileName.cpp, HeaderName.h
 - Class/structure name: ClassName
 - Variable name: variableName
 - Function/method name: this_is_a_function_name

   AMReX API exception (including Grid and GridAccess):
   - Methods inherited from AmrCore/AmrMesh keep their original AMReX names,
     including Geom(), boxArray(), DistributionMap(), refRatio(),
     finestLevel(), maxLevel(), gridEff(), and SetGridEff().
   - Grid uses these inherited methods directly. Do not redefine them or add
     snake_case aliases merely to adapt their spelling to FLEKS conventions.
   - GridAccess forwards needed AMReX queries with the same original names,
     for example Geom() and gridEff(). It remains a read-only query facade;
     this naming exception does not authorize adding mutation methods.
   - Required AMReX virtual overrides retain their exact library names and
     signatures. Calls on other AMReX objects also keep library names.
   - Methods introduced by FLEKS follow snake_case, including Grid-owned
     queries such as node_box_array() and output methods such as write_mf().
     Existing helpers that add behavior, such as n_lev(), are FLEKS methods,
     not mere spelling aliases for inherited methods.
   - Preserve overloads, constness, return semantics, and ownership when
     applying naming changes. Identify the owning API before renaming a call.

3. 'using namespace amrex' is allowed in *.cpp files. Otherwise, do NOT leave 
   'using namespace xxx' in the code. 

4. For beauty, the order of headers: std headers -> AMReX headers -> user headers. 
    The behavior of the code shoud NOT depend on the order of headers.      

5. Use 'nullptr' instead of NULL for pointer. 

6. Always use 'const' if possible. 

7. Lambda is useful, but a regular function is better if it is universal or it is long. 

8. Using debug flags and Valgrind to check errors.

9. Follow the coventional commits format: https://www.conventionalcommits.org/en/v1.0.0/
