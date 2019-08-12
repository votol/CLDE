__kernel void Daxpy(__global double* y, __global double* x, const double coe) {

   int gid = get_global_id(0);
   
   y[gid] += coe * x[gid];
}
