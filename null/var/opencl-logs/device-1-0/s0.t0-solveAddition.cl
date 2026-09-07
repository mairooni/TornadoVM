#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void solveAddition(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *matrix1, __global uchar *matrix2, __global uchar *result)
{
  long l_26, l_27; 
  double d_29, d_12, d_7, d_37, d_33; 
  ulong ul_11, ul_16, ul_24, ul_23, ul_30, ul_28, ul_2, ul_34, ul_1, ul_0, ul_32, ul_31, ul_6, ul_5, ul_4, ul_36, ul_3, ul_35, ul_10, ul_9, ul_8; 
  int i_14, i_13, i_25, i_39, i_22, i_38, i_21, i_20, i_19, i_18, i_17, i_15; 

  // BLOCK 0
  ul_0  =  (ulong) matrix1;
  ul_1  =  (ulong) matrix2;
  ul_2  =  (ulong) result;
  ul_3  =  ul_0 + 32L;
  ul_4  =  *((__global ulong *) ul_3);
  ul_5  =  ul_0 + ul_4;
  ul_6  =  ul_5 + 24L;
  d_7  =  *((__global double *) ul_6);
  ul_8  =  ul_1 + 32L;
  ul_9  =  *((__global ulong *) ul_8);
  ul_10  =  ul_1 + ul_9;
  ul_11  =  ul_10 + 96L;
  d_12  =  *((__global double *) ul_11);
  i_13  =  isequal(d_7, d_12);
  if(i_13 == 1)
  {
    // BLOCK 1
    i_14  =  get_global_size(0);
    i_15  =  get_global_size(1);
    ul_16  =  ul_2 + 32L;
    i_17  =  get_global_id(0);
    i_18  =  get_global_id(1);
    // BLOCK 2 MERGES [1 6 ]
    i_19  =  i_18;
    for(;i_19 < 8;)
    {
      // BLOCK 3
      i_20  =  i_19 << 3;
      i_21  =  i_20 + 3;
      // BLOCK 4 MERGES [3 5 ]
      i_22  =  i_17;
      for(;i_22 < 8;)
      {
        // BLOCK 5
        ul_23  =  *((__global ulong *) ul_3);
        ul_24  =  ul_0 + ul_23;
        i_25  =  i_21 + i_22;
        l_26  =  (long) i_25;
        l_27  =  l_26 << 3;
        ul_28  =  ul_24 + l_27;
        d_29  =  *((__global double *) ul_28);
        ul_30  =  *((__global ulong *) ul_8);
        ul_31  =  ul_1 + ul_30;
        ul_32  =  ul_31 + l_27;
        d_33  =  *((__global double *) ul_32);
        ul_34  =  *((__global ulong *) ul_16);
        ul_35  =  ul_2 + ul_34;
        ul_36  =  ul_35 + l_27;
        d_37  =  d_29 + d_33;
        *((__global double *) ul_36)  =  d_37;
        i_38  =  i_14 + i_22;
        i_22  =  i_38;
      }  // B5
      // BLOCK 6
      i_39  =  i_15 + i_19;
      i_19  =  i_39;
    }  // B6
    // BLOCK 7
    return;
    else
    {
      // BLOCK 8
      return;
    }  // B8
  }  //  kernel
