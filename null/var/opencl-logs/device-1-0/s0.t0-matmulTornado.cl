#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void matmulTornado(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *thisx, __global uchar *that, __global uchar *out, __private int dim1)
{
  int i_45, i_46, i_43, i_50, i_29, i_28, i_31, i_30, i_35, i_13, i_12, i_15, i_14, i_9, i_73, i_8, i_72, i_11, i_10, i_19, i_63, i_57, i_59, i_58, i_4, i_68, i_7, i_6, i_3; 
  long l_65, l_16, l_64, l_20, l_74, l_75, l_60; 
  ulong ul_21, ul_2, ul_66, ul_1, ul_17, ul_0, ul_61, ul_76; 
  bool b_44, b_33, b_32, b_37, b_41, b_55, b_39, b_38; 
  short sh_27; 
  uint ui_26, ui_25, ui_24, ui_23; 
  char ch_18, ch_62, ch_22; 
  float f_53, f_54, f_56, f_49, f_51, f_52, f_47, f_48, f_42, f_5, f_69, f_70, f_71, f_40, f_34, f_67, f_36; 

  // BLOCK 0
  ul_0  =  (ulong) thisx;
  ul_1  =  (ulong) that;
  ul_2  =  (ulong) out;
  i_3  =  get_global_id(0);
  i_4  =  i_3 << 13;
  // BLOCK 1 MERGES [0 20 ]
  f_5  =  0.0F;
  i_6  =  0;
  for(;i_6 < 8192;)
  {
    // BLOCK 2
    i_7  =  i_4 + i_6;
    i_8  =  i_7 >> 31;
    i_9  =  i_8 >> 27;
    i_10  =  i_9 + i_7;
    i_11  =  i_10 >> 5;
    i_12  =  i_11 << 1;
    i_13  =  i_10 & -32;
    i_14  =  i_12 + i_13;
    i_15  =  i_14 + 24;
    l_16  =  (long) i_15;
    ul_17  =  ul_0 + l_16;
    ch_18  =  *((__global char *) ul_17);
    i_19  =  i_14 + 25;
    l_20  =  (long) i_19;
    ul_21  =  ul_0 + l_20;
    ch_22  =  *((__global char *) ul_21);
    ui_23  =  ch_22 & 255U;
    ui_24  =  ui_23 << 8;
    ui_25  =  ch_18 & 255U;
    ui_26  =  ui_24 | ui_25;
    sh_27  =  (short) ui_26;
    i_28  =  (int) sh_27;
    i_29  =  i_28 & 32768;
    i_30  =  i_28 & 31744;
    i_31  =  i_30 >> 10;
    b_32  =  i_31 == 31;
    if(b_32)
    {
      // BLOCK 3
      b_33  =  (i_29 & -32768) == 0;
      if(b_33)
      {
        // BLOCK 4
        f_34  =  InfinityF;
      }  // B4
      else
      {
        // BLOCK 5
        f_34  =  -InfinityF;
      }  // B5
    }  // B3
    else
    {
      // BLOCK 6
      i_35  =  i_28 & 1023;
      f_36  =  (float) i_35;
      b_37  =  (i_30 & -1024) == 0;
      if(b_37)
      {
        // BLOCK 7
        b_38  =  (i_28 & 1023) == 0;
        if(b_38)
        {
          // BLOCK 8
          b_39  =  (i_29 & -32768) == 0;
          if(b_39)
          {
            // BLOCK 9
            f_34  =  0.0F;
          }  // B9
          else
          {
            // BLOCK 10
            f_34  =  -0.0F;
          }  // B10
        }  // B8
        else
        {
          // BLOCK 11
          f_40  =  f_36 * 5.9604645E-8F;
          b_41  =  (i_29 & -32768) == 0;
          if(b_41)
          {
            // BLOCK 12
            f_34  =  f_40;
          }  // B12
          else
          {
            // BLOCK 13
            f_42  =  -f_40;
            f_34  =  f_42;
          }  // B13
        }  // B11
      }  // B7
      else
      {
        // BLOCK 14
        i_43  =  i_31 + -15;
        b_44  =  i_31 < 15;
        if(b_44)
        {
          // BLOCK 15
          i_45  =  -i_43;
          i_46  =  1 << i_45;
          f_47  =  (float) i_46;
          f_48  =  1.0F / f_47;
          f_49  =  f_48;
        }  // B15
        else
        {
          // BLOCK 16
          i_50  =  1 << i_43;
          f_51  =  (float) i_50;
          f_49  =  f_51;
        }  // B16
        // BLOCK 17 MERGES [16 15 ]
        f_52  =  f_36 / 1024.0F;
        f_53  =  f_52 + 1.0F;
        f_54  =  f_53 * f_49;
        b_55  =  (i_29 & -32768) == 0;
        if(b_55)
        {
          // BLOCK 18
          f_34  =  f_54;
        }  // B18
        else
        {
          // BLOCK 19
          f_56  =  -f_54;
          f_34  =  f_56;
        }  // B19
      }  // B17
    }  // B6
    // BLOCK 20 MERGES [4 9 12 18 5 10 13 19 ]
    i_57  =  i_7 % 32;
    i_58  =  i_57 + i_14;
    i_59  =  i_58 + 26;
    l_60  =  (long) i_59;
    ul_61  =  ul_0 + l_60;
    ch_62  =  *((__global char *) ul_61);
    i_63  =  i_6 + 6;
    l_64  =  (long) i_63;
    l_65  =  l_64 << 2;
    ul_66  =  ul_1 + l_65;
    f_67  =  *((__global float *) ul_66);
    i_68  =  (int) ch_62;
    f_69  =  (float) i_68;
    f_70  =  f_34 * f_69;
    f_71  =  fma(f_70, f_67, f_5);
    i_72  =  i_6 + 1;
    f_5  =  f_71;
    i_6  =  i_72;
  }  // B20
  // BLOCK 21
  i_73  =  i_3 + 6;
  l_74  =  (long) i_73;
  l_75  =  l_74 << 2;
  ul_76  =  ul_2 + l_75;
  *((__global float *) ul_76)  =  f_5;
  return;
}  //  kernel
