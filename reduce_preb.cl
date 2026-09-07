#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void reductionOneBlockWithLayer(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *output, __global uchar *x, __private int size, __private float ermsNorm, __private int localMemSize)
{
  int i_6, i_4, i_3, i_29, i_28, i_25, i_22, i_19, i_18, i_17, i_16, i_15, i_14; 
  long l_30, l_31, l_7, l_8; 
  float f_12, f_13, f_10, f_11, f_21, f_24, f_23, f_27, f_37, f_35, f_41, f_39, f_45, f_43, f_49, f_47, f_52, f_53, f_51, f_56, f_57, f_54, f_55, f_60, f_61, f_58, f_59; 
  bool b_33, b_5, b_20, b_26; 
  ulong ul_0, ul_32, ul_1, ul_34, ul_36, ul_38, ul_40, ul_9, ul_42, ul_44, ul_46, ul_48, ul_50; 

  // BLOCK 0
  ul_0  =  (ulong) output;
  ul_1  =  (ulong) x;
  __local float adf_2[128];
  i_3  =  get_local_id(0);
  i_4  =  get_global_id(0);
  b_5  =  i_4 < 1024;
  if(b_5)
  {
    // BLOCK 1
    i_6  =  i_4 + 6;
    l_7  =  (long) i_6;
    l_8  =  l_7 << 2;
    ul_9  =  ul_1 + l_8;
    f_10  =  *((__global float *) ul_9);
    adf_2[i_3]  =  f_10;
    f_11  =  adf_2[i_3];
    f_12  =  adf_2[i_3];
    f_13  =  f_11 * f_12;
    adf_2[i_3]  =  f_13;
  }  // B1
  else
  {
    // BLOCK 2
    adf_2[i_3]  =  0.0F;
  }  // B2
  barrier(CLK_LOCAL_MEM_FENCE);
  printf(">> barrier\n");
  // BLOCK 3 MERGES [1 2 ]
  i_14  =  get_local_size(0);
  i_15  =  i_14 >> 31;
  i_16  =  i_15 + i_14;
  i_17  =  i_16 >> 1;
  // BLOCK 4 MERGES [3 8 ]
  i_18  =  i_17;
  for(;i_18 >= 1;)
  {
    // BLOCK 5
    barrier(CLK_LOCAL_MEM_FENCE);
    i_19  =  i_18 >> 1;
    b_20  =  i_3 < i_18;
    if(b_20)
    {
      // BLOCK 6
      f_21  =  adf_2[i_3];
      i_22  =  i_18 + i_3;
      f_23  =  adf_2[i_22];
      f_24  =  f_21 + f_23;
      adf_2[i_3]  =  f_24;
    }  // B6
    else
    {
      // BLOCK 7
    }  // B7
    // BLOCK 8 MERGES [7 6 ]
    i_25  =  i_19;
    i_18  =  i_25;
  }  // B8
  // BLOCK 9
  b_26  =  i_3 == 0;
  if(b_26)
  {
    // BLOCK 10
    f_27  =  adf_2[0];
    i_28  =  get_group_id(0);
    i_29  =  i_28 + 7;
    l_30  =  (long) i_29;
    l_31  =  l_30 << 2;
    ul_32  =  ul_0 + l_31;
    *((__global float *) ul_32)  =  f_27;
  }  // B10
  else
  {
    // BLOCK 11
  }  // B11
  // BLOCK 12 MERGES [11 10 ]
  b_33  =  i_4 == 0;
  if(b_33)
  {
    // BLOCK 13
    ul_34  =  ul_0 + 28L;
    f_35  =  *((__global float *) ul_34);
    ul_36  =  ul_0 + 32L;
    f_37  =  *((__global float *) ul_36);
    ul_38  =  ul_0 + 36L;
    f_39  =  *((__global float *) ul_38);
    ul_40  =  ul_0 + 40L;
    f_41  =  *((__global float *) ul_40);
    ul_42  =  ul_0 + 44L;
    f_43  =  *((__global float *) ul_42);
    ul_44  =  ul_0 + 48L;
    f_45  =  *((__global float *) ul_44);
    ul_46  =  ul_0 + 52L;
    f_47  =  *((__global float *) ul_46);
    ul_48  =  ul_0 + 56L;
    f_49  =  *((__global float *) ul_48);
    ul_50  =  ul_0 + 24L;
    f_51  =  f_35 + 0.0F;
    f_52  =  f_51 + f_37;
    f_53  =  f_52 + f_39;
    f_54  =  f_53 + f_41;
    f_55  =  f_54 + f_43;
    f_56  =  f_55 + f_45;
    f_57  =  f_56 + f_47;
    f_58  =  f_57 + f_49;
    f_59  =  f_58 / 1024.0F;
    f_60  =  f_59 + 1.0E-5F;
    f_61  =  rsqrt(f_60);
    *((__global float *) ul_50)  =  f_61;
  }  // B13
  else
  {
    // BLOCK 14
  }  // B14
  // BLOCK 15 MERGES [14 13 ]
  return;
}  //  kernel

