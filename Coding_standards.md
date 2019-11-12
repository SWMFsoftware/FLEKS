## Visit the God of C++: https://isocpp.org/wiki/faq

## Standards
1. Avoid using raw pointer and 'new' to allocate memory unless there is a good reason. 
 Try to use shared_ptr or unique_ptr to manage resources. For arrays, using an 
 object to manage resources is a good choice. 

2. Naming standards:
 - File name: FileName.cpp, HeaderName.h
 - Class/structure name: ClassName
 - Variable name: variableName
 - Function/method name: this_is_a_function_name

3. Never put 'using namespace xxx' into a header. 

4. For beauty, the order of headers: std headers -> AMReX headers -> user headers. 
    The correctness of the code shoud NOT depend on the order of headers.      

4. Use 'nullptr' instead of NULL for pointer. 

6. Always use 'const' if possible. 

7. Lambda is useful, but a regular function is better if it is universal or it is long. 

8. Using debug flags and Valgrind to check errors. 
