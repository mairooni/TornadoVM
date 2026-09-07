#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void computeInit(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *this)
{
  long l_11, l_10; 
  int i_9, i_15, i_1, i_6, i_5; 
  ulong ul_14, ul_0, ul_3, ul_4, ul_7, ul_8, ul_13; 
  uint ui_12, ui_2; 

  // BLOCK 0
  ul_0  =  (ulong) this;
  i_1  =  get_global_size(0);
  ul_4  =  ul_0 + 12L;
  ui_2  =  *((__global uint *) ul_4);
  ul_3  =  ul_0 ((ulong) ui_2 << 3);
  i_5  =  get_global_id(0);
  // BLOCK 1 MERGES [0 2 ]
  i_6  =  i_5;
  for(;i_6 < 1024;)
  {
    // BLOCK 2
    ul_7  =  *((__global ulong *) ul_3);
    ul_8  =  ul_0 + ul_7;
    i_9  =  i_6 + 4;
    l_10  =  (long) i_9;
    l_11  =  l_10 << 2;
    ul_14  =  ul_8 + l_11;
    ui_12  =  *((__global uint *) ul_14);
    ul_13  =  ul_8 ((ulong) ui_12 << 3);
    *((__global int *) ul_13)  =  100;
    i_15  =  i_1 + i_6;
    i_6  =  i_15;
  }  // B2
  // BLOCK 3
  return;
}  //  kernel
