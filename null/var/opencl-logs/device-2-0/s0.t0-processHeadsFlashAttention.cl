#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void processHeadsFlashAttention(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *q, __global uchar *key_cache, __global uchar *value_cache, __global uchar *xb, __private int nHeads, __private int headSize, __private int kvDim, __private int kvMul, __global uchar *positionHolder, __private int layer, __private int contextLength)
{
  long l_24, l_25, l_50, l_114, l_49, l_113; 
  int i_12, i_18, i_17, i_16, i_15, i_22, i_21, i_20, i_19, i_106, i_105, i_103, i_110, i_107, i_112, i_117, i_90, i_87, i_93, i_98, i_97, i_100, i_74, i_72, i_78, i_75, i_82, i_84, i_58, i_57, i_56, i_62, i_60, i_59, i_64, i_70, i_69, i_67, i_42, i_41, i_40, i_39, i_46, i_45, i_44, i_43, i_48, i_47, i_53, i_23, i_30, i_28, i_34, i_33, i_32, i_31, i_38, i_35; 
  ulong ul_14, ul_0, ul_1, ul_2, ul_3, ul_51, ul_115, ul_4, ul_54, ul_26; 
  bool b_13, b_79; 
  float f_109, f_111, f_108, f_37, f_101, f_102, f_104, f_99, f_36, f_61, f_63, f_55, f_52, f_116, f_77, f_80, f_73, f_76, f_71, f_65, f_66, f_68, f_29, f_94, f_95, f_96, f_89, f_27, f_91, f_92, f_85, f_86, f_88, f_81, f_83; 

  // BLOCK 0
  ul_0  =  (ulong) q;
  ul_1  =  (ulong) key_cache;
  ul_2  =  (ulong) value_cache;
  ul_3  =  (ulong) xb;
  ul_4  =  (ulong) positionHolder;
  __private float ul_5[512];
  __private float* ul_6 = ul_5;
  __local float adf_7[512];
  __local float adf_8[2048];
  __local float adf_9[2048];
  __local float adf_10[4];
  __local float adf_11[1];
  i_12  =  get_group_id(0);
  b_13  =  i_12 < 32;
  if(b_13)
  {
    // BLOCK 1
    ul_14  =  ul_4 + 24L;
    i_15  =  *((__global int *) ul_14);
    // BLOCK 2 MERGES [1 3 ]
    i_16  =  0;
    for(;i_16 < 512;)
    {
      // BLOCK 3
      ul_6[i_16]  =  0.0F;
      i_17  =  i_16 + 1;
      i_16  =  i_17;
    }  // B3
    // BLOCK 4
    i_18  =  get_local_size(0);
    i_19  =  i_12 << 9;
    i_20  =  i_19 + 6;
    i_21  =  get_local_id(0);
    // BLOCK 5 MERGES [4 6 ]
    i_22  =  i_21;
    for(;i_22 < 512;)
    {
      // BLOCK 6
      i_23  =  i_20 + i_22;
      l_24  =  (long) i_23;
      l_25  =  l_24 << 2;
      ul_26  =  ul_0 + l_25;
      f_27  =  *((__global float *) ul_26);
      adf_7[i_22]  =  f_27;
      i_28  =  i_18 + i_22;
      i_22  =  i_28;
    }  // B6
    // BLOCK 7
    barrier(CLK_LOCAL_MEM_FENCE);
    f_29  =  -1.0F / 0.0F;
    i_30  =  i_12 >> 31;
    i_31  =  i_30 >> 30;
    i_32  =  i_31 + i_12;
    i_33  =  i_32 >> 2;
    i_34  =  i_33 << 9;
    i_35  =  i_34 + 6;
    // BLOCK 8 MERGES [7 44 ]
    f_36  =  f_29;
    f_37  =  0.0F;
    i_38  =  0;
    for(;i_15 >= i_38;)
    {
      // BLOCK 9
      i_39  =  i_38 + 3;
      i_40  =  min(i_39, i_15);
      i_41  =  i_38 + i_21;
      // BLOCK 10 MERGES [9 14 ]
      i_42  =  i_41;
      for(;i_40 >= i_42;)
      {
        // BLOCK 11
        i_43  =  i_42 - i_38;
        i_44  =  i_43 << 9;
        i_45  =  i_42 << 12;
        i_46  =  i_45 + i_35;
        // BLOCK 12 MERGES [11 13 ]
        i_47  =  0;
        for(;i_47 < 512;)
        {
          // BLOCK 13
          i_48  =  i_46 + i_47;
          l_49  =  (long) i_48;
          l_50  =  l_49 << 2;
          ul_51  =  ul_1 + l_50;
          f_52  =  *((__global float *) ul_51);
          i_53  =  i_44 + i_47;
          adf_8[i_53]  =  f_52;
          ul_54  =  ul_2 + l_50;
          f_55  =  *((__global float *) ul_54);
          adf_9[i_53]  =  f_55;
          i_56  =  i_47 + 1;
          i_47  =  i_56;
        }  // B13
        // BLOCK 14
        i_57  =  i_42 + i_18;
        i_42  =  i_57;
      }  // B14
      // BLOCK 15
      barrier(CLK_LOCAL_MEM_FENCE);
      // BLOCK 16 MERGES [15 20 ]
      i_58  =  i_41;
      for(;i_40 >= i_58;)
      {
        // BLOCK 17
        i_59  =  i_58 - i_38;
        i_60  =  i_59 << 9;
        // BLOCK 18 MERGES [17 19 ]
        f_61  =  0.0F;
        i_62  =  0;
        for(;i_62 < 512;)
        {
          // BLOCK 19
          f_63  =  adf_7[i_62];
          i_64  =  i_60 + i_62;
          f_65  =  adf_8[i_64];
          f_66  =  fma(f_63, f_65, f_61);
          i_67  =  i_62 + 1;
          f_61  =  f_66;
          i_62  =  i_67;
        }  // B19
        // BLOCK 20
        f_68  =  f_61 / 22.627417F;
        adf_10[i_59]  =  f_68;
        i_69  =  i_58 + i_18;
        i_58  =  i_69;
      }  // B20
      // BLOCK 21
      barrier(CLK_LOCAL_MEM_FENCE);
      i_70  =  i_40 - i_38;
      // BLOCK 22 MERGES [21 26 ]
      f_71  =  f_29;
      i_72  =  0;
      for(;i_70 >= i_72;)
      {
        // BLOCK 23
        f_73  =  adf_10[i_72];
        i_74  =  i_72 + 1;
        i_75  =  isless(f_71, f_73);
        if(i_75 == 1)
        {
          // BLOCK 24
          f_76  =  adf_10[i_72];
          f_77  =  f_76;
        }  // B24
        else
        {
          // BLOCK 25
          f_77  =  f_71;
        }  // B25
        // BLOCK 26 MERGES [25 24 ]
        i_78  =  i_74;
        f_71  =  f_77;
        i_72  =  i_78;
      }  // B26
      // BLOCK 27
      b_79  =  i_21 == 0;
      if(b_79)
      {
        // BLOCK 28
        adf_11[0]  =  f_71;
      }  // B28
      else
      {
        // BLOCK 29
      }  // B29
      // BLOCK 30 MERGES [29 28 ]
      barrier(CLK_LOCAL_MEM_FENCE);
      f_80  =  adf_11[0];
      f_81  =  fmax(f_36, f_80);
      i_82  =  isequal(f_36, f_81);
      if(i_82 == 1)
      {
        // BLOCK 31
        f_83  =  f_37;
      }  // B31
      else
      {
        // BLOCK 32
        i_84  =  isequal(f_36, f_29);
        if(i_84 == 1)
        {
          // BLOCK 33
          f_83  =  f_37;
        }  // B33
        else
        {
          // BLOCK 34
          f_85  =  f_36 - f_81;
          f_86  =  exp(f_85);
          // BLOCK 35 MERGES [34 36 ]
          i_87  =  0;
          for(;i_87 < 512;)
          {
            // BLOCK 36
            f_88  =  ul_6[i_87];
            f_89  =  f_86 * f_88;
            ul_6[i_87]  =  f_89;
            i_90  =  i_87 + 1;
            i_87  =  i_90;
          }  // B36
          // BLOCK 37
          f_91  =  f_37 * f_86;
          f_83  =  f_91;
        }  // B34
      }  // B32
      // BLOCK 38 MERGES [31 33 37 ]
      f_92  =  f_83;
      // BLOCK 39 MERGES [38 43 ]
      f_92  =  f_83;
      i_93  =  0;
      for(;i_70 >= i_93;)
      {
        // BLOCK 40
        f_94  =  adf_10[i_93];
        f_95  =  f_94 - f_81;
        f_96  =  exp(f_95);
        i_97  =  i_93 << 9;
        // BLOCK 41 MERGES [40 42 ]
        i_98  =  0;
        for(;i_98 < 512;)
        {
          // BLOCK 42
          f_99  =  ul_6[i_98];
          i_100  =  i_97 + i_98;
          f_101  =  adf_9[i_100];
          f_102  =  fma(f_96, f_101, f_99);
          ul_6[i_98]  =  f_102;
          i_103  =  i_98 + 1;
          i_98  =  i_103;
        }  // B42
        // BLOCK 43
        f_104  =  f_92 + f_96;
        i_105  =  i_93 + 1;
        f_92  =  f_104;
        i_93  =  i_105;
      }  // B43
      // BLOCK 44
      barrier(CLK_LOCAL_MEM_FENCE);
      i_106  =  i_38 + 4;
      f_36  =  f_81;
      f_37  =  f_92;
      i_38  =  i_106;
      // BLOCK 45
      i_107  =  isless(0.0F, f_37);
      if(i_107 == 1)
      {
        // BLOCK 46
        f_108  =  1.0F / f_37;
        f_109  =  f_108;
      }  // B46
      else
      {
        // BLOCK 47
        f_109  =  0.0F;
      }  // B47
      // BLOCK 48 MERGES [46 47 ]
      // BLOCK 49 MERGES [48 50 ]
      i_110  =  i_21;
      for(;i_110 < 512;)
      {
        // BLOCK 50
        f_111  =  ul_6[i_110];
        i_112  =  i_110 + i_20;
        l_113  =  (long) i_112;
        l_114  =  l_113 << 2;
        ul_115  =  ul_3 + l_114;
        f_116  =  f_109 * f_111;
        *((__global float *) ul_115)  =  f_116;
        i_117  =  i_110 + i_18;
        i_110  =  i_117;
      }  // B50
      // BLOCK 51
      return;
    }  // B1
    else
    {
      // BLOCK 52
      return;
    }  // B52
  }  //  kernel
