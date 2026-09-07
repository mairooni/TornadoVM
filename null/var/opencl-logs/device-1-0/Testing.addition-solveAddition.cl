#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void solveAddition(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *matrix1, __global uchar *matrix2, __global uchar *result)
{
  long l_34, l_33; 
  double d_15, d_17, d_36, d_21, d_7, d_40, d_12, d_44; 
  int i_22, i_24, i_23, i_26, i_25, i_28, i_27, i_14, i_46, i_29, i_45, i_32; 
  ulong ul_6, ul_4, ul_5, ul_10, ul_11, ul_8, ul_9, ul_2, ul_3, ul_0, ul_1, ul_20, ul_13, ul_18, ul_19, ul_16, ul_38, ul_39, ul_37, ul_42, ul_43, ul_41, ul_30, ul_31, ul_35; 

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
  ul_13  =  ul_2 + 32L;
  i_14  =  isequal(d_7, d_12);
  if(i_14 == 1)
  {
    // BLOCK 1
  }  // B1
  else
  {
    // BLOCK 2
    d_15  =  *((__global double *) ul_6);
    ul_16  =  ul_10 + 24L;
    d_17  =  *((__global double *) ul_16);
    ul_18  =  *((__global ulong *) ul_13);
    ul_19  =  ul_2 + ul_18;
    ul_20  =  ul_19 + 24L;
    d_21  =  d_15 + d_17;
    *((__global double *) ul_20)  =  d_21;
  }  // B2
  // BLOCK 3 MERGES [1 2 ]
  i_22  =  get_global_size(0);
  i_23  =  get_global_size(1);
  i_24  =  get_global_id(0);
  i_25  =  get_global_id(1);
  // BLOCK 4 MERGES [3 8 ]
  i_26  =  i_25;
  for(;i_26 < 8;)
  {
    // BLOCK 5
    i_27  =  i_26 << 3;
    i_28  =  i_27 + 3;
    // BLOCK 6 MERGES [5 7 ]
    i_29  =  i_24;
    for(;i_29 < 8;)
    {
      // BLOCK 7
      ul_30  =  *((__global ulong *) ul_3);
      ul_31  =  ul_0 + ul_30;
      i_32  =  i_28 + i_29;
      l_33  =  (long) i_32;
      l_34  =  l_33 << 3;
      ul_35  =  ul_31 + l_34;
      d_36  =  *((__global double *) ul_35);
      ul_37  =  *((__global ulong *) ul_8);
      ul_38  =  ul_1 + ul_37;
      ul_39  =  ul_38 + l_34;
      d_40  =  *((__global double *) ul_39);
      ul_41  =  *((__global ulong *) ul_13);
      ul_42  =  ul_2 + ul_41;
      ul_43  =  ul_42 + l_34;
      d_44  =  d_36 + d_40;
      *((__global double *) ul_43)  =  d_44;
      i_45  =  i_22 + i_29;
      i_29  =  i_45;
    }  // B7
    // BLOCK 8
    i_46  =  i_23 + i_26;
    i_26  =  i_46;
  }  // B8
  // BLOCK 9
  return;
}  // B9
}  //  kernel
