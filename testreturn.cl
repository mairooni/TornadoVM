#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void conditionAfterInnerForLoopReturn(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *matrixA, __global uchar *matrixB, __global uchar *result)
{
  int i_11, i_12, i_28, i_15, i_3, i_4, i_5, i_37, i_6, i_38, i_7, i_39; 
  ulong ul_26, ul_25, ul_24, ul_31, ul_30, ul_29, ul_18, ul_22, ul_21, ul_20, ul_10, ul_9, ul_8, ul_14, ul_13, ul_35, ul_2, ul_34, ul_1, ul_33, ul_0; 
  long l_16, l_17; 
  double d_32, d_27, d_23, d_36, d_19; 

  // BLOCK 0
  ul_0  =  (ulong) matrixA;
  ul_1  =  (ulong) matrixB;
  ul_2  =  (ulong) result;
  i_3  =  get_global_id(1);
  // BLOCK 1 MERGES [0 6 ]
  i_4  =  i_3;
  for(;i_4 < 8;)
  {
    // BLOCK 2
    i_5  =  i_4 << 3;
    i_6  =  i_5 + 3;
    i_7  =  get_global_size(0);
    ul_8  =  ul_2 + 32L;
    ul_9  =  ul_1 + 32L;
    ul_10  =  ul_0 + 32L;
    i_11  =  get_global_id(0);
    // BLOCK 3 MERGES [2 4 ]
    i_12  =  i_11;
    for(;i_12 < 8;)
    {
      // BLOCK 4
      ul_13  =  *((__global ulong *) ul_10);
      ul_14  =  ul_0 + ul_13;
      i_15  =  i_6 + i_12;
      l_16  =  (long) i_15;
      l_17  =  l_16 << 3;
      ul_18  =  ul_14 + l_17;
      d_19  =  *((__global double *) ul_18);
      ul_20  =  *((__global ulong *) ul_9);
      ul_21  =  ul_1 + ul_20;
      ul_22  =  ul_21 + l_17;
      d_23  =  *((__global double *) ul_22);
      ul_24  =  *((__global ulong *) ul_8);
      ul_25  =  ul_2 + ul_24;
      ul_26  =  ul_25 + l_17;
      d_27  =  d_19 + d_23;
      *((__global double *) ul_26)  =  d_27;
      i_28  =  i_7 + i_12;
      i_12  =  i_28;
    }  // B4
    // BLOCK 5
    ul_29  =  *((__global ulong *) ul_10);
    ul_30  =  ul_0 + ul_29;
    ul_31  =  ul_30 + 24L;
    d_32  =  *((__global double *) ul_31);
    ul_33  =  *((__global ulong *) ul_9);
    ul_34  =  ul_1 + ul_33;
    ul_35  =  ul_34 + 24L;
    d_36  =  *((__global double *) ul_35);
    i_37  =  isequal(d_32, d_36);
    if(i_37 == 1)
    {
      // BLOCK 6
      printf(">> d_32: %lf is equal to d_36: %lf, continue\n", d_32, d_36);
      i_38  =  get_global_size(1);
      i_39  =  i_38 + i_4;
      i_4  =  i_39;
    }  // B6
    else
    {
      // BLOCK 7
      printf("$$ d_32: %lf is NOT equal to d_36: %lf, return\n", d_32, d_36);
      return;
      printf("This should never print if return works\n");
    }  // B7
    // BLOCK 8
    return;
  }  // B8
}  //  kernel
