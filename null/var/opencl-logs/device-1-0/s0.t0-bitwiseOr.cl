#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void bitwiseOr(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *result, __global uchar *input, __global uchar *elements)
{
  char ch_5, ch_19, ch_35, ch_16, ch_31, ch_12, ch_10, ch_24, ch_25; 
  long l_29, l_22, l_8; 
  int i_27, i_26, i_21, i_20, i_34, i_33, i_32, i_28, i_7, i_6, i_4, i_36, i_18, i_17, i_15, i_14, i_13; 
  ulong ul_23, ul_11, ul_9, ul_0, ul_30, ul_3, ul_1, ul_2; 

  // BLOCK 0
  ul_0  =  (ulong) result;
  ul_1  =  (ulong) input;
  ul_2  =  (ulong) elements;
  ul_3  =  ul_2 + 24L;
  // BLOCK 1 MERGES [0 2 ]
  i_4  =  0;
  for(;;)
  {
    ch_5  =  *((__global char *) ul_3);
    i_6  =  (int) ch_5;
    if(!(i_4 < i_6))
    {
      break;
    }
    // BLOCK 2
    i_7  =  i_4 + 24;
    l_8  =  (long) i_7;
    ul_9  =  ul_0 + l_8;
    ch_10  =  *((__global char *) ul_9);
    ul_11  =  ul_1 + l_8;
    ch_12  =  *((__global char *) ul_11);
    i_13  =  (int) ch_10;
    i_14  =  (int) ch_12;
    i_15  =  i_13 | i_14;
    ch_16  =  (char) i_15;
    *((__global char *) ul_9)  =  ch_16;
    i_17  =  i_4 + 1;
    i_4  =  i_17;
  }  // B2
  // BLOCK 3
  // BLOCK 4 MERGES [3 5 ]
  i_18  =  0;
  for(;;)
  {
    ch_19  =  *((__global char *) ul_3);
    i_20  =  (int) ch_19;
    if(!(i_18 < i_20))
    {
      break;
    }
    // BLOCK 5
    i_21  =  i_18 + 24;
    l_22  =  (long) i_21;
    ul_23  =  ul_0 + l_22;
    ch_24  =  *((__global char *) ul_23);
    ch_25  =  *((__global char *) ul_3);
    i_26  =  (int) ch_25;
    i_27  =  i_26 + i_18;
    i_28  =  i_27 + 24;
    l_29  =  (long) i_28;
    ul_30  =  ul_1 + l_29;
    ch_31  =  *((__global char *) ul_30);
    i_32  =  (int) ch_24;
    i_33  =  (int) ch_31;
    i_34  =  i_32 | i_33;
    ch_35  =  (char) i_34;
    *((__global char *) ul_23)  =  ch_35;
    i_36  =  i_18 + 1;
    i_18  =  i_36;
  }  // B5
  // BLOCK 6
  return;
}  // B3
}  //  kernel
