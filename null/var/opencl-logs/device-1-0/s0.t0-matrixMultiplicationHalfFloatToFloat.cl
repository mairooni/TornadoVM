#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void matrixMultiplicationHalfFloatToFloat(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *A, __global uchar *B, __global uchar *C, __private int size)
{
  int i_20, i_14, i_15, i_28, i_21, i_34, i_3, i_35, i_4, i_29, i_9, i_10, i_11, i_12, i_5, i_6, i_7, i_8; 
  ulong ul_1, ul_0, ul_32, ul_2, ul_18, ul_24; 
  float f_33; 
  long l_23, l_22, l_17, l_31, l_16, l_30; 
  half half_19, half_25, half_13, half_26, half_27; 

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
    i_9  =  i_8 + 6;
    i_10  =  i_8 + 12;
    // BLOCK 3 MERGES [2 7 ]
    i_11  =  i_5;
    for(;i_11 < 256;)
    {
      // BLOCK 4
      i_12  =  i_11 + 12;
      // BLOCK 5 MERGES [4 6 ]
      half_13  =  0.0;
      i_14  =  0;
      for(;i_14 < 256;)
      {
        // BLOCK 6
        i_15  =  i_10 + i_14;
        l_16  =  (long) i_15;
        l_17  =  l_16 << 1;
        ul_18  =  ul_0 + l_17;
        half_19  =  *((__global half *) ul_18);
        i_20  =  i_14 << 8;
        i_21  =  i_20 + i_12;
        l_22  =  (long) i_21;
        l_23  =  l_22 << 1;
        ul_24  =  ul_1 + l_23;
        half_25  =  *((__global half *) ul_24);
        half_26  =  half_19 * half_25;
        half_27  =  half_13 + half_26;
        i_28  =  i_14 + 1;
        half_13  =  half_27;
        i_14  =  i_28;
      }  // B6
      // BLOCK 7
      i_29  =  i_9 + i_11;
      l_30  =  (long) i_29;
      l_31  =  l_30 << 2;
      ul_32  =  ul_2 + l_31;
      f_33  =  convert_float(half_13);
      *((__global float *) ul_32)  =  f_33;
      i_34  =  i_3 + i_11;
      i_11  =  i_34;
    }  // B7
    // BLOCK 8
    i_35  =  i_4 + i_7;
    i_7  =  i_35;
  }  // B8
  // BLOCK 9
  return;
}  //  kernel
