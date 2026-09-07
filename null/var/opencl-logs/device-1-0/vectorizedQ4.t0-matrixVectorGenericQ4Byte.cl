#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
float uk_ac_manchester_tornado_examples_compute_MatrixVectorRowMajor_matrixVectorRowMajorOptimizedQ4_0Byte(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, ulong context, int localSize, ulong x, ulong q, int n)
{
  float f_94, f_95, f_96, f_97, f_161, f_90, f_154, f_91, f_92, f_93, f_102, f_166, f_103, f_104, f_98, f_99, f_163, f_100, f_164, f_110, f_72, f_14, f_15, f_16, f_13, f_77, f_150, f_87, f_151, f_88, f_152, f_89, f_153, f_82; 
  int i_73, i_68, i_67, i_66, i_64, i_63, i_62, i_61, i_83, i_78, i_40, i_39, i_38, i_37, i_165, i_36, i_35, i_162, i_159, i_30, i_158, i_29, i_157, i_28, i_156, i_59, i_58, i_57, i_56, i_54, i_53, i_52, i_51, i_49, i_48, i_44, i_11, i_139, i_10, i_138, i_8, i_7, i_6, i_134, i_5, i_133, i_4, i_132, i_3, i_131, i_2, i_130, i_1, i_129, i_124, i_27, i_155, i_26, i_25, i_24, i_23, i_22, i_21, i_20, i_19, i_18, i_146, i_17, i_145, i_144, i_142, i_141, i_12, i_140, i_107, i_106, i_105, i_101, i_123, i_122, i_121, i_120, i_119, i_118, i_117, i_116, i_115, i_114, i_113, i_112, i_111, i_109, i_108; 
  long l_85, l_84, l_148, l_147, l_80, l_79, l_45, l_75, l_74, l_41, l_135, l_70, l_69, l_32, l_31, l_126, l_125; 
  bool b_65, b_50, b_55, b_60, b_143, b_160; 
  ulong ul_33, ul_81, ul_71, ul_86, ul_149, ul_42, ul_9, ul_136, ul_127, ul_46, ul_76; 
  half half_128, half_34; 
  char ch_47, ch_137, ch_43; 

  // BLOCK 0
  __local float adf_0[localSize];
  i_1  =  localSize << 2;
  i_2  =  n + 31;
  i_3  =  i_2 >> 31;
  i_4  =  i_3 >> 27;
  i_5  =  i_4 + i_2;
  i_6  =  i_5 >> 5;
  i_7  =  get_group_id(0);
  i_8  =  i_6 * i_7;
  ul_9  =  q + 24L;
  i_10  =  get_local_id(0);
  i_11  =  i_10 << 2;
  i_12  =  n + -3;
  // BLOCK 1 MERGES [0 14 ]
  f_13  =  0.0F;
  f_14  =  0.0F;
  f_15  =  0.0F;
  f_16  =  0.0F;
  i_17  =  i_11;
  for(;i_17 < i_12;)
  {
    // BLOCK 2
    i_18  =  *((__global int *) ul_9);
    i_19  =  i_17 >> 31;
    i_20  =  i_19 >> 27;
    i_21  =  i_20 + i_17;
    i_22  =  i_21 >> 5;
    i_23  =  i_8 + i_22;
    i_24  =  i_23 << 4;
    i_25  =  i_23 << 1;
    i_26  =  i_24 + i_25;
    i_27  =  i_26 + i_18;
    i_28  =  i_27 >> 31;
    i_29  =  i_28 + i_27;
    i_30  =  i_29 >> 1;
    l_31  =  (long) i_30;
    l_32  =  l_31 << 1;
    ul_33  =  q + l_32;
    half_34  =  *((__global half *) ul_33);
    i_35  =  i_17 % 32;
    i_36  =  i_35 >> 31;
    i_37  =  i_36 + i_35;
    i_38  =  i_37 >> 1;
    i_39  =  i_38 + i_26;
    i_40  =  i_39 + 18;
    l_41  =  (long) i_40;
    ul_42  =  q + l_41;
    ch_43  =  *((__global char *) ul_42);
    i_44  =  i_39 + 19;
    l_45  =  (long) i_44;
    ul_46  =  q + l_45;
    ch_47  =  *((__global char *) ul_46);
    i_48  =  (int) ch_43;
    i_49  =  i_48 & 15;
    b_50  =  i_49 >= 0 && i_49 < 8;
    if(b_50)
    {
      // BLOCK 3
      i_51  =  i_49;
    }  // B3
    else
    {
      // BLOCK 4
      i_52  =  i_49 + -16;
      i_51  =  i_52;
    }  // B4
    // BLOCK 5 MERGES [3 4 ]
    i_53  =  i_48 >> 4;
    i_54  =  i_53 & 15;
    b_55  =  i_54 >= 0 && i_54 < 8;
    if(b_55)
    {
      // BLOCK 6
      i_56  =  i_54;
    }  // B6
    else
    {
      // BLOCK 7
      i_57  =  i_54 + -16;
      i_56  =  i_57;
    }  // B7
    // BLOCK 8 MERGES [6 7 ]
    i_58  =  (int) ch_47;
    i_59  =  i_58 & 15;
    b_60  =  i_59 >= 0 && i_59 < 8;
    if(b_60)
    {
      // BLOCK 9
      i_61  =  i_59;
    }  // B9
    else
    {
      // BLOCK 10
      i_62  =  i_59 + -16;
      i_61  =  i_62;
    }  // B10
    // BLOCK 11 MERGES [9 10 ]
    i_63  =  i_58 >> 4;
    i_64  =  i_63 & 15;
    b_65  =  i_64 >= 0 && i_64 < 8;
    if(b_65)
    {
      // BLOCK 12
      i_66  =  i_64;
    }  // B12
    else
    {
      // BLOCK 13
      i_67  =  i_64 + -16;
      i_66  =  i_67;
    }  // B13
    // BLOCK 14 MERGES [12 13 ]
    i_68  =  i_17 + 4;
    l_69  =  (long) i_68;
    l_70  =  l_69 << 2;
    ul_71  =  x + l_70;
    f_72  =  *((__global float *) ul_71);
    i_73  =  i_17 + 5;
    l_74  =  (long) i_73;
    l_75  =  l_74 << 2;
    ul_76  =  x + l_75;
    f_77  =  *((__global float *) ul_76);
    i_78  =  i_17 + 6;
    l_79  =  (long) i_78;
    l_80  =  l_79 << 2;
    ul_81  =  x + l_80;
    f_82  =  *((__global float *) ul_81);
    i_83  =  i_17 + 7;
    l_84  =  (long) i_83;
    l_85  =  l_84 << 2;
    ul_86  =  x + l_85;
    f_87  =  *((__global float *) ul_86);
    f_88  =  convert_float((float) half_34);
    f_89  =  (float) i_66;
    f_90  =  f_88 * f_89;
    f_91  =  fma(f_90, f_87, f_16);
    f_92  =  (float) i_61;
    f_93  =  f_92 * f_88;
    f_94  =  fma(f_93, f_82, f_15);
    f_95  =  (float) i_56;
    f_96  =  f_95 * f_88;
    f_97  =  fma(f_96, f_77, f_14);
    f_98  =  (float) i_51;
    f_99  =  f_98 * f_88;
    f_100  =  fma(f_99, f_72, f_13);
    i_101  =  i_1 + i_17;
    f_13  =  f_100;
    f_14  =  f_97;
    f_15  =  f_94;
    f_16  =  f_91;
    i_17  =  i_101;
  }  // B14
  // BLOCK 15
  f_102  =  f_13 + f_14;
  f_103  =  f_102 + f_15;
  f_104  =  f_103 + f_16;
  i_105  =  n >> 31;
  i_106  =  i_105 >> 30;
  i_107  =  i_106 + n;
  i_108  =  i_107 & -4;
  i_109  =  i_108 + i_10;
  // BLOCK 16 MERGES [15 20 ]
  f_110  =  f_104;
  i_111  =  i_109;
  for(;i_111 < n;)
  {
    // BLOCK 17
    i_112  =  *((__global int *) ul_9);
    i_113  =  i_111 >> 31;
    i_114  =  i_113 >> 27;
    i_115  =  i_114 + i_111;
    i_116  =  i_115 >> 5;
    i_117  =  i_116 + i_8;
    i_118  =  i_117 << 4;
    i_119  =  i_117 << 1;
    i_120  =  i_118 + i_119;
    i_121  =  i_120 + i_112;
    i_122  =  i_121 >> 31;
    i_123  =  i_122 + i_121;
    i_124  =  i_123 >> 1;
    l_125  =  (long) i_124;
    l_126  =  l_125 << 1;
    ul_127  =  q + l_126;
    half_128  =  *((__global half *) ul_127);
    i_129  =  i_111 % 32;
    i_130  =  i_129 >> 31;
    i_131  =  i_130 + i_129;
    i_132  =  i_131 >> 1;
    i_133  =  i_132 + i_120;
    i_134  =  i_133 + 18;
    l_135  =  (long) i_134;
    ul_136  =  q + l_135;
    ch_137  =  *((__global char *) ul_136);
    i_138  =  (int) ch_137;
    i_139  =  i_129 & 1;
    i_140  =  i_139 << 2;
    i_141  =  i_138 >> i_140;
    i_142  =  i_141 & 15;
    b_143  =  i_142 >= 0 && i_142 < 8;
    if(b_143)
    {
      // BLOCK 18
      i_144  =  i_142;
    }  // B18
    else
    {
      // BLOCK 19
      i_145  =  i_142 + -16;
      i_144  =  i_145;
    }  // B19
    // BLOCK 20 MERGES [18 19 ]
    i_146  =  i_111 + 4;
    l_147  =  (long) i_146;
    l_148  =  l_147 << 2;
    ul_149  =  x + l_148;
    f_150  =  *((__global float *) ul_149);
    f_151  =  convert_float((float) half_128);
    f_152  =  (float) i_144;
    f_153  =  f_151 * f_152;
    f_154  =  fma(f_153, f_150, f_110);
    i_155  =  localSize + i_111;
    f_110  =  f_154;
    i_111  =  i_155;
  }  // B20
  // BLOCK 21
  adf_0[i_10]  =  f_110;
  barrier(CLK_LOCAL_MEM_FENCE);
  i_156  =  localSize >> 31;
  i_157  =  i_156 + localSize;
  i_158  =  i_157 >> 1;
  // BLOCK 22 MERGES [21 26 ]
  i_159  =  i_158;
  for(;i_159 >= 1;)
  {
    // BLOCK 23
    b_160  =  i_10 < i_159;
    if(b_160)
    {
      // BLOCK 24
      f_161  =  adf_0[i_10];
      i_162  =  i_159 + i_10;
      f_163  =  adf_0[i_162];
      f_164  =  f_161 + f_163;
      adf_0[i_10]  =  f_164;
    }  // B24
    else
    {
      // BLOCK 25
    }  // B25
    // BLOCK 26 MERGES [25 24 ]
    barrier(CLK_LOCAL_MEM_FENCE);
    i_165  =  i_159 >> 1;
    i_159  =  i_165;
  }  // B26
  // BLOCK 27
  f_166  =  adf_0[0];
  return f_166;
}  //  kernel

#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void matrixVectorGenericQ4Byte(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *x, __global uchar *output, __global uchar *q, __private int dim1, __private int dim0, __private int localWorkGroupSize)
{
  float f_7; 
  int i_9, i_6, i_4; 
  bool b_8, b_5; 
  long l_11, l_10; 
  ulong ul_3, ul_2, ul_1, ul_0, ul_12; 

  // BLOCK 0
  ul_0  =  (ulong) context;
  ul_1  =  (ulong) x;
  ul_2  =  (ulong) output;
  ul_3  =  (ulong) q;
  i_4  =  get_group_id(0);
  b_5  =  i_4 < 2048;
  if(b_5)
  {
    // BLOCK 1
    i_6  =  get_local_id(0);
    f_7  =  uk_ac_manchester_tornado_examples_compute_MatrixVectorRowMajor_matrixVectorRowMajorOptimizedQ4_0Byte(_kernel_context, _constant_region, _local_region, _atomics, ul_0, 32, ul_1, ul_3, 8192);
    b_8  =  i_6 == 0;
    if(b_8)
    {
      // BLOCK 2
      i_9  =  i_4 + 4;
      l_10  =  (long) i_9;
      l_11  =  l_10 << 2;
      ul_12  =  ul_2 + l_11;
      *((__global float *) ul_12)  =  f_7;
      return;
    }  // B2
    else
    {
      // BLOCK 3
      return;
    }  // B3
    // BLOCK 4 MERGES [3 2 ]
    return;
    else
    {
      // BLOCK 5
      return;
    }  // B5
  }  //  kernel
