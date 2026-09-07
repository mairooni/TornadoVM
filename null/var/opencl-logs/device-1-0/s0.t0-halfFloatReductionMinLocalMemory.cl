#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void halfFloatReductionMinLocalMemory(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *a, __global uchar *b)
{
  long l_29, l_12, l_30, l_11, l_10, l_5, l_6; 
  ulong ul_1, ul_0, ul_31, ul_7; 
  half half_20, half_22, half_23, half_8, half_26; 
  int i_28, i_27, i_24, i_21, i_18, i_16, i_17, i_14, i_15, i_13, i_9, i_4, i_3; 
  bool b_19, b_25; 

  // BLOCK 0
  ul_0  =  (ulong) a;
  ul_1  =  (ulong) b;
  __local half half_2[256];
  i_3  =  get_global_id(0);
  i_4  =  i_3 + 8;
  l_5  =  (long) i_4;
  l_6  =  l_5 << 1;
  ul_7  =  ul_0 + l_6;
  half_8  =  *((__global half *) ul_7);
  i_9  =  get_local_id(0);
  l_10  =  (long) i_9;
  l_11  =  l_10 << 2;
  l_12  =  l_11 + 16L;
  half_2[i_9]  =  half_8;
  i_13  =  get_local_size(0);
  i_14  =  i_13 >> 31;
  i_15  =  i_14 + i_13;
  i_16  =  i_15 >> 1;
  // BLOCK 1 MERGES [0 5 ]
  i_17  =  i_16;
  for(;i_17 >= 1;)
  {
    // BLOCK 2
    barrier(CLK_LOCAL_MEM_FENCE);
    i_18  =  i_17 >> 1;
    b_19  =  i_9 < i_17;
    if(b_19)
    {
      // BLOCK 3
      half_20  =  half_2[i_9];
      i_21  =  i_17 + i_9;
      half_22  =  half_2[i_21];
      half_23  =  fmin(half_20, half_22);
      half_2[i_9]  =  half_23;
    }  // B3
    else
    {
      // BLOCK 4
    }  // B4
    // BLOCK 5 MERGES [4 3 ]
    i_24  =  i_18;
    i_17  =  i_24;
  }  // B5
  // BLOCK 6
  b_25  =  i_9 == 0;
  if(b_25)
  {
    // BLOCK 7
    half_26  =  half_2[0];
    i_27  =  get_group_id(0);
    i_28  =  i_27 + 8;
    l_29  =  (long) i_28;
    l_30  =  l_29 << 1;
    ul_31  =  ul_1 + l_30;
    *((__global half *) ul_31)  =  half_26;
  }  // B7
  else
  {
    // BLOCK 8
  }  // B8
  // BLOCK 9 MERGES [8 7 ]
  return;
}  //  kernel
