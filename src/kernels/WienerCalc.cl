__kernel void wienerCalc(__global double* out, __global double* random, const unsigned int shift, const double time) 
{
    __local double tmp[32]; 
    int id = get_global_id(0);
    int loc_id = get_local_id(0);
    int gid = id / 32;
    int ind_a = (shift + loc_id) % 97;
    int ind_b = (64 + shift + loc_id) % 97;
    tmp[loc_id] = random[97 * gid + ind_a] - random[97 * gid + ind_b];
    if(tmp[loc_id] < 0.0)
        tmp[loc_id] += 1.0;
    random[97*gid + ind_a] = tmp[loc_id];
    if(id < DIM)
    {
        double u1 = tmp[loc_id - loc_id%2] - 0.25*convert_double(loc_id%2);
        double u2 = tmp[loc_id - loc_id%2 + 1];
        out[id] = sqrt(-2.0 * time *  log(u2)) * cos(2.0 * 3.14159265358979323846 * u1);
    }
}
