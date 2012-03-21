

__kernel void gravlens_kernel(
                     __global float* input,
                     __global float* output,
                     const unsigned int count
                    )
{
   unsigned int i = get_global_id(0);
   if(i < count)
       output[i] = input[i] - input[i];
}

