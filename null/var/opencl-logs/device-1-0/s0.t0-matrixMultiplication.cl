#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void matrixMultiplication(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *A, __global uchar *B, __global uchar *C, __private int size)
{
  half half_18, half_24, half_25; 
  long l_29, l_30, l_15, l_16, l_21, l_22; 
  ulong ul_23, ul_2, ul_0, ul_1, ul_17, ul_31; 
  int i_20, i_26, i_27, i_12, i_13, i_14, i_19, i_4, i_5, i_6, i_7, i_8, i_9, i_10, i_11, i_28, i_32, i_33, i_3; 

  // BLOCK 0
  ul_0  =  (ulong) A;
  ul_1  =  (ulong) B;
  ul_2  =  (ulong) C;
  i_3  =  get_global_size(0);
  i_4  =  get_global_size(1);
  i_5  =  get_global_id(0);
  i_6  =  get_global_id(1);
  // BLOCK 1 MERGES [0 8 ]
  i_7  =  i_6;
  for(;i_7 < 256;)
  {
    // BLOCK 2
    i_8  =  i_7 << 8;
    i_9  =  i_8 + 12;
    // BLOCK 3 MERGES [2 7 ]
    i_10  =  i_5;
    for(;i_10 < 256;)
    {
      // BLOCK 4
      i_11  =  i_10 + 12;
      // BLOCK 5 MERGES [4 6 ]
      i_12  =  0.0;
      i_13  =  0;
      for(;i_13 < 256;)
      {
        // BLOCK 6
        i_14  =  i_9 + i_13;
        l_15  =  (long) i_14;
        l_16  =  l_15 << 1;
        ul_17  =  ul_0 + l_16;
        half_18  =  *((__global half *) ul_17);
        i_19  =  i_13 << 8;
        i_20  =  i_19 + i_11;
        l_21  =  (long) i_20;
        l_22  =  l_21 << 1;
        ul_23  =  ul_1 + l_22;
        half_24  =  *((__global half *) ul_23);
        half_25  =  half_18 * half_24;
        i_26  =  i_12 + half_25;
        i_27  =  i_13 + 1;
        i_12  =  i_26;
        i_13  =  i_27;
      }  // B6
      // BLOCK 7
      i_28  =  i_10 + i_9;
      l_29  =  (long) i_28;
      l_30  =  l_29 << 1;
      ul_31  =  ul_2 + l_30;
      *((__global half *) ul_31)  =  i_12;
      i_32  =  i_3 + i_10;
      i_10  =  i_32;
    }  // B7
    // BLOCK 8
    i_33  =  i_4 + i_7;
    i_7  =  i_33;
  }  // B8
  // BLOCK 9
  return;
}  //  kernel
