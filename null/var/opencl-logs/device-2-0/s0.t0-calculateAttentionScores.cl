#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void calculateAttentionScores(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *positionNlayer, __private int seqLen, __global uchar *query, __global uchar *keyCache, __global uchar *attScores, __private int kvDim, __private int kvMul, __private int headSize, __private int loff, __private int localWorkgroupSize)
{
  ulong ul_31, ul_0, ul_1, ul_2, ul_3, ul_4, ul_38, ul_26; 
  long l_25, l_30, l_29, l_36, l_37, l_24; 
  float f_33, f_32, f_27, f_21, f_39; 
  int i_7, i_8, i_40, i_5, i_6, i_11, i_12, i_9, i_10, i_35, i_34, i_23, i_22, i_28, i_15, i_16, i_13, i_14, i_19, i_20, i_17, i_18; 

  // BLOCK 0
  ul_0  =  (ulong) positionNlayer;
  ul_1  =  (ulong) query;
  ul_2  =  (ulong) keyCache;
  ul_3  =  (ulong) attScores;
  ul_4  =  ul_0 + 24L;
  i_5  =  *((__global int *) ul_4);
  i_6  =  get_local_size(0);
  i_7  =  get_group_id(0);
  i_8  =  i_7 << 1;
  i_9  =  i_8 + 6;
  i_10  =  i_7 / 6;
  i_11  =  i_10 << 3;
  i_12  =  i_11 - i_10;
  i_13  =  i_12 + 9;
  i_14  =  i_7 << 3;
  i_15  =  i_14 - i_7;
  i_16  =  i_15 + 6;
  i_17  =  get_local_id(0);
  // BLOCK 1 MERGES [0 5 ]
  i_18  =  i_17;
  for(;i_5 >= i_18;)
  {
    // BLOCK 6
    return;
  }  // B6
  // BLOCK 2
  i_19  =  i_18 << 2;
  i_20  =  i_19 + i_13;
  // BLOCK 3 MERGES [2 4 ]
  f_21  =  0.0F;
  i_22  =  0;
  for(;i_22 < 8192;)
  {
    // BLOCK 4
    i_23  =  i_16 + i_22;
    l_24  =  (long) i_23;
    l_25  =  l_24 << 2;
    ul_26  =  ul_1 + l_25;
    f_27  =  *((__global float *) ul_26);
    i_28  =  i_20 + i_22;
    l_29  =  (long) i_28;
    l_30  =  l_29 << 2;
    ul_31  =  ul_2 + l_30;
    f_32  =  *((__global float *) ul_31);
    f_33  =  fma(f_27, f_32, f_21);
    i_34  =  i_22 + 1;
    f_21  =  f_33;
    i_22  =  i_34;
  }  // B4
  // BLOCK 5
  i_35  =  i_9 + i_18;
  l_36  =  (long) i_35;
  l_37  =  l_36 << 2;
  ul_38  =  ul_3 + l_37;
  f_39  =  f_21 / 2.6457512F;
  *((__global float *) ul_38)  =  f_39;
  i_40  =  i_6 + i_18;
  i_18  =  i_40;
}  // B5
}  //  kernel
