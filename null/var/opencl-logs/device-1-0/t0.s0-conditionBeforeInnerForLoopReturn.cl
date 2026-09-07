#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void conditionBeforeInnerForLoopReturn(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *matrixA, __global uchar *matrixB, __global uchar *result)
{
  int i_25, i_15, i_16, i_17, i_18, i_3, i_19, i_4, i_21, i_22, i_38, i_39; 
  double d_37, d_33, d_29, d_14, d_9; 
  long l_27, l_26; 
  ulong ul_5, ul_6, ul_7, ul_8, ul_10, ul_11, ul_12, ul_30, ul_31, ul_0, ul_32, ul_1, ul_2, ul_34, ul_35, ul_36, ul_23, ul_24, ul_28, ul_13, ul_20; 

  // BLOCK 0
  ul_0  =  (ulong) matrixA;
  ul_1  =  (ulong) matrixB;
  ul_2  =  (ulong) result;
  i_3  =  get_global_id(1);
  // BLOCK 1 MERGES [0 6 ]
  i_4  =  i_3;
  for(;i_4 < 128;)
  {
    // BLOCK 2
    ul_5  =  ul_0 + 32L;
    ul_6  =  *((__global ulong *) ul_5);
    ul_7  =  ul_0 + ul_6;
    ul_8  =  ul_7 + 24L;
    d_9  =  *((__global double *) ul_8);
    ul_10  =  ul_1 + 32L;
    ul_11  =  *((__global ulong *) ul_10);
    ul_12  =  ul_1 + ul_11;
    ul_13  =  ul_12 + 1056L;
    d_14  =  *((__global double *) ul_13);
    i_15  =  isequal(d_9, d_14);
    if(i_15 == 1)
    {
      // BLOCK 3
      i_16  =  i_4 << 7;
      i_17  =  i_16 + 3;
      i_18  =  get_global_size(0);
      i_19  =  get_global_size(1);
      ul_20  =  ul_2 + 32L;
      i_21  =  get_global_id(0);
      // BLOCK 4 MERGES [3 5 ]
      i_22  =  i_21;
      for(;i_22 < 128;)
      {
        // BLOCK 5
        ul_23  =  *((__global ulong *) ul_5);
        ul_24  =  ul_0 + ul_23;
        i_25  =  i_17 + i_22;
        l_26  =  (long) i_25;
        l_27  =  l_26 << 3;
        ul_28  =  ul_24 + l_27;
        d_29  =  *((__global double *) ul_28);
        ul_30  =  *((__global ulong *) ul_10);
        ul_31  =  ul_1 + ul_30;
        ul_32  =  ul_31 + l_27;
        d_33  =  *((__global double *) ul_32);
        ul_34  =  *((__global ulong *) ul_20);
        ul_35  =  ul_2 + ul_34;
        ul_36  =  ul_35 + l_27;
        d_37  =  d_29 + d_33;
        *((__global double *) ul_36)  =  d_37;
        i_38  =  i_18 + i_22;
        i_22  =  i_38;
      }  // B5
      // BLOCK 6
      i_39  =  i_19 + i_4;
      i_4  =  i_39;
    }  // B6
    else
    {
      // BLOCK 7
      return;
    }  // B7
    // BLOCK 8
    return;
  }  //  kernel
