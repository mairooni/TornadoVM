#pragma OPENCL EXTENSION cl_khr_fp64 : enable  
#pragma OPENCL EXTENSION cl_khr_fp16 : enable  
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics : enable  
__kernel void precompute_faulty_calculation_data_serial(__global long *_kernel_context, __constant uchar *_constant_region, __local uchar *_local_region, __global int *_atomics, __global uchar *initialCoords, __global uchar *precomputedArgsAndUpperTri, __private int N)
{
  long l_192, l_254, l_253, l_116, l_244, l_115, l_243, l_177, l_55, l_54, l_76, l_204, l_138, l_266, l_137, l_265, l_144, l_272, l_143, l_271, l_77, l_205, l_132, l_260, l_131, l_259, l_193, l_71, l_199, l_70, l_198, l_89, l_83, l_211, l_82, l_210, l_88, l_150, l_149, l_42, l_41, l_176, l_47, l_46; 
  ulong ul_141, ul_147, ul_146, ul_145, ul_152, ul_151, ul_154, ul_153, ul_168, ul_167, ul_166, ul_171, ul_170, ul_169, ul_190, ul_189, ul_196, ul_195, ul_194, ul_200, ul_202, ul_201, ul_208, ul_207, ul_206, ul_212, ul_215, ul_214, ul_213, ul_228, ul_227, ul_232, ul_231, ul_230, ul_229, ul_251, ul_250, ul_0, ul_256, ul_255, ul_4, ul_1, ul_257, ul_7, ul_263, ul_6, ul_262, ul_5, ul_261, ul_268, ul_11, ul_267, ul_10, ul_9, ul_15, ul_14, ul_13, ul_269, ul_276, ul_19, ul_275, ul_18, ul_274, ul_17, ul_273, ul_23, ul_22, ul_21, ul_27, ul_26, ul_25, ul_288, ul_31, ul_30, ul_29, ul_35, ul_34, ul_290, ul_33, ul_289, ul_37, ul_68, ul_67, ul_72, ul_74, ul_73, ul_80, ul_79, ul_78, ul_84, ul_86, ul_85, ul_92, ul_91, ul_90, ul_93, ul_108, ul_107, ul_106, ul_105, ul_110, ul_109, ul_128, ul_129, ul_135, ul_134, ul_133, ul_140, ul_139; 
  int i_81, i_209, i_148, i_136, i_264, i_75, i_203, i_142, i_270, i_291, i_38, i_87, i_112, i_240, i_111, i_175, i_239, i_50, i_114, i_242, i_113, i_241, i_52, i_51, i_53, i_40, i_39, i_172, i_174, i_238, i_45, i_173, i_191, i_130, i_258, i_69, i_197, i_252; 
  double d_246, d_119, d_247, d_120, d_248, d_121, d_249, d_122, d_124, d_125, d_126, d_127, d_99, d_100, d_101, d_102, d_103, d_104, d_234, d_235, d_236, d_237, d_216, d_217, d_218, d_219, d_220, d_221, d_94, d_222, d_95, d_223, d_96, d_224, d_97, d_225, d_98, d_226, d_180, d_181, d_182, d_183, d_185, d_58, d_186, d_59, d_187, d_60, d_188, d_61, d_63, d_64, d_65, d_66, d_163, d_36, d_164, d_165, d_43, d_44, d_48, d_49, d_20, d_277, d_278, d_279, d_24, d_280, d_281, d_282, d_155, d_283, d_28, d_156, d_284, d_157, d_285, d_158, d_286, d_159, d_287, d_32, d_160, d_161, d_162, d_8, d_12, d_16; 
  bool b_117, b_245, b_178, b_179, b_62, b_123, b_56, b_184, b_57, b_233, b_118; 

  // BLOCK 0
  ul_0  =  (ulong) initialCoords;
  ul_1  =  (ulong) precomputedArgsAndUpperTri;
  __private double ul_2[8];
  __private double* ul_3 = ul_2;
  ul_4  =  ul_0 + 32L;
  ul_5  =  *((__global ulong *) ul_4);
  ul_6  =  ul_0 + ul_5;
  ul_7  =  ul_6 + 24L;
  d_8  =  *((__global double *) ul_7);
  ul_3[24L]  =  d_8;
  ul_9  =  *((__global ulong *) ul_4);
  ul_10  =  ul_0 + ul_9;
  ul_11  =  ul_10 + 32L;
  d_12  =  *((__global double *) ul_11);
  ul_3[32L]  =  d_12;
  ul_13  =  *((__global ulong *) ul_4);
  ul_14  =  ul_0 + ul_13;
  ul_15  =  ul_14 + 40L;
  d_16  =  *((__global double *) ul_15);
  ul_3[40L]  =  d_16;
  ul_17  =  *((__global ulong *) ul_4);
  ul_18  =  ul_0 + ul_17;
  ul_19  =  ul_18 + 48L;
  d_20  =  *((__global double *) ul_19);
  ul_3[48L]  =  d_20;
  ul_21  =  *((__global ulong *) ul_4);
  ul_22  =  ul_0 + ul_21;
  ul_23  =  ul_22 + 56L;
  d_24  =  *((__global double *) ul_23);
  ul_3[56L]  =  d_24;
  ul_25  =  *((__global ulong *) ul_4);
  ul_26  =  ul_0 + ul_25;
  ul_27  =  ul_26 + 64L;
  d_28  =  *((__global double *) ul_27);
  ul_3[64L]  =  d_28;
  ul_29  =  *((__global ulong *) ul_4);
  ul_30  =  ul_0 + ul_29;
  ul_31  =  ul_30 + 72L;
  d_32  =  *((__global double *) ul_31);
  ul_3[72L]  =  d_32;
  ul_33  =  *((__global ulong *) ul_4);
  ul_34  =  ul_0 + ul_33;
  ul_35  =  ul_34 + 80L;
  d_36  =  *((__global double *) ul_35);
  ul_3[80L]  =  d_36;
  ul_37  =  ul_1 + 32L;
  // BLOCK 1 MERGES [0 32 ]
  i_38  =  0;
  for(;i_38 < 4;)
  {
    // BLOCK 2
    i_39  =  i_38 << 1;
    i_40  =  i_39 + 3;
    l_41  =  (long) i_40;
    l_42  =  l_41 << 3;
    d_43  =  ul_3[l_42];
    d_44  =  d_43 / 57.29577951308232;
    ul_3[l_42]  =  d_44;
    i_45  =  i_39 + 4;
    l_46  =  (long) i_45;
    l_47  =  l_46 << 3;
    d_48  =  ul_3[l_47];
    d_49  =  d_48 / 57.29577951308232;
    ul_3[l_47]  =  d_49;
    i_50  =  i_38 << 2;
    i_51  =  i_38 << 4;
    i_52  =  i_50 + i_51;
    i_53  =  i_52 + 7;
    l_54  =  (long) i_53;
    l_55  =  l_54 << 3;
    b_56  =  i_38 < 1;
    if(b_56)
    {
      // BLOCK 3
      b_57  =  i_38 == 0;
      if(b_57)
      {
        // BLOCK 4
        d_58  =  ul_3[24L];
        d_59  =  d_58 / 57.29577951308232;
        ul_3[24L]  =  d_59;
        d_60  =  ul_3[32L];
        d_61  =  d_60 / 57.29577951308232;
        ul_3[32L]  =  d_61;
      }  // B4
      else
      {
        // BLOCK 5
      }  // B5
      // BLOCK 6 MERGES [5 4 ]
      b_62  =  i_38 < 0;
      if(b_62)
      {
        // BLOCK 7
        d_63  =  ul_3[l_42];
        d_64  =  ul_3[l_47];
        d_65  =  ul_3[24L];
        d_66  =  ul_3[32L];
        ul_67  =  *((__global ulong *) ul_37);
        ul_68  =  ul_1 + ul_67;
        i_69  =  i_52 + 3;
        l_70  =  (long) i_69;
        l_71  =  l_70 << 3;
        ul_72  =  ul_68 + l_71;
        *((__global double *) ul_72)  =  d_63;
        ul_73  =  *((__global ulong *) ul_37);
        ul_74  =  ul_1 + ul_73;
        i_75  =  i_52 + 4;
        l_76  =  (long) i_75;
        l_77  =  l_76 << 3;
        ul_78  =  ul_74 + l_77;
        *((__global double *) ul_78)  =  d_64;
        ul_79  =  *((__global ulong *) ul_37);
        ul_80  =  ul_1 + ul_79;
        i_81  =  i_52 + 5;
        l_82  =  (long) i_81;
        l_83  =  l_82 << 3;
        ul_84  =  ul_80 + l_83;
        *((__global double *) ul_84)  =  d_65;
        ul_85  =  *((__global ulong *) ul_37);
        ul_86  =  ul_1 + ul_85;
        i_87  =  i_52 + 6;
        l_88  =  (long) i_87;
        l_89  =  l_88 << 3;
        ul_90  =  ul_86 + l_89;
        *((__global double *) ul_90)  =  d_66;
        ul_91  =  *((__global ulong *) ul_37);
        ul_92  =  ul_1 + ul_91;
        ul_93  =  ul_92 + l_55;
        d_94  =  native_sin(d_63);
        d_95  =  native_sin(d_65);
        d_96  =  native_cos(d_63);
        d_97  =  native_cos(d_65);
        d_98  =  d_96 * d_97;
        d_99  =  d_64 - d_66;
        d_100  =  native_cos(d_99);
        d_101  =  d_98 * d_100;
        d_102  =  fma(d_94, d_95, d_101);
        d_103  =  acos(d_102);
        d_104  =  d_103 * 6366707.019493707;
        *((__global double *) ul_93)  =  d_104;
      }  // B7
      else
      {
        // BLOCK 8
        ul_105  =  *((__global ulong *) ul_37);
        ul_106  =  ul_1 + ul_105;
        ul_107  =  ul_106 + l_55;
        *((__global double *) ul_107)  =  0.0;
      }  // B8
    }  // B6
    else
    {
      // BLOCK 9
      ul_108  =  *((__global ulong *) ul_37);
      ul_109  =  ul_1 + ul_108;
      ul_110  =  ul_109 + l_55;
      *((__global double *) ul_110)  =  0.0;
    }  // B9
    // BLOCK 10 MERGES [7 9 8 ]
    i_111  =  i_50 + 1;
    i_112  =  i_111 << 2;
    i_113  =  i_112 + i_111;
    i_114  =  i_113 + 7;
    l_115  =  (long) i_114;
    l_116  =  l_115 << 3;
    b_117  =  i_38 < 2;
    if(b_117)
    {
      // BLOCK 11
      b_118  =  i_38 == 0;
      if(b_118)
      {
        // BLOCK 12
        d_119  =  ul_3[40L];
        d_120  =  d_119 / 57.29577951308232;
        ul_3[40L]  =  d_120;
        d_121  =  ul_3[48L];
        d_122  =  d_121 / 57.29577951308232;
        ul_3[48L]  =  d_122;
      }  // B12
      else
      {
        // BLOCK 13
      }  // B13
      // BLOCK 14 MERGES [13 12 ]
      b_123  =  i_38 < 1;
      if(b_123)
      {
        // BLOCK 15
        d_124  =  ul_3[l_42];
        d_125  =  ul_3[l_47];
        d_126  =  ul_3[40L];
        d_127  =  ul_3[48L];
        ul_128  =  *((__global ulong *) ul_37);
        ul_129  =  ul_1 + ul_128;
        i_130  =  i_113 + 3;
        l_131  =  (long) i_130;
        l_132  =  l_131 << 3;
        ul_133  =  ul_129 + l_132;
        *((__global double *) ul_133)  =  d_124;
        ul_134  =  *((__global ulong *) ul_37);
        ul_135  =  ul_1 + ul_134;
        i_136  =  i_113 + 4;
        l_137  =  (long) i_136;
        l_138  =  l_137 << 3;
        ul_139  =  ul_135 + l_138;
        *((__global double *) ul_139)  =  d_125;
        ul_140  =  *((__global ulong *) ul_37);
        ul_141  =  ul_1 + ul_140;
        i_142  =  i_113 + 5;
        l_143  =  (long) i_142;
        l_144  =  l_143 << 3;
        ul_145  =  ul_141 + l_144;
        *((__global double *) ul_145)  =  d_126;
        ul_146  =  *((__global ulong *) ul_37);
        ul_147  =  ul_1 + ul_146;
        i_148  =  i_113 + 6;
        l_149  =  (long) i_148;
        l_150  =  l_149 << 3;
        ul_151  =  ul_147 + l_150;
        *((__global double *) ul_151)  =  d_127;
        ul_152  =  *((__global ulong *) ul_37);
        ul_153  =  ul_1 + ul_152;
        ul_154  =  ul_153 + l_116;
        d_155  =  native_sin(d_124);
        d_156  =  native_sin(d_126);
        d_157  =  native_cos(d_124);
        d_158  =  native_cos(d_126);
        d_159  =  d_157 * d_158;
        d_160  =  d_125 - d_127;
        d_161  =  native_cos(d_160);
        d_162  =  d_159 * d_161;
        d_163  =  fma(d_155, d_156, d_162);
        d_164  =  acos(d_163);
        d_165  =  d_164 * 6366707.019493707;
        *((__global double *) ul_154)  =  d_165;
      }  // B15
      else
      {
        // BLOCK 16
        ul_166  =  *((__global ulong *) ul_37);
        ul_167  =  ul_1 + ul_166;
        ul_168  =  ul_167 + l_116;
        *((__global double *) ul_168)  =  0.0;
      }  // B16
    }  // B14
    else
    {
      // BLOCK 17
      ul_169  =  *((__global ulong *) ul_37);
      ul_170  =  ul_1 + ul_169;
      ul_171  =  ul_170 + l_116;
      *((__global double *) ul_171)  =  0.0;
    }  // B17
    // BLOCK 18 MERGES [15 17 16 ]
    i_172  =  i_50 + 2;
    i_173  =  i_172 << 2;
    i_174  =  i_173 + i_172;
    i_175  =  i_174 + 7;
    l_176  =  (long) i_175;
    l_177  =  l_176 << 3;
    b_178  =  i_38 < 3;
    if(b_178)
    {
      // BLOCK 19
      b_179  =  i_38 == 0;
      if(b_179)
      {
        // BLOCK 20
        d_180  =  ul_3[56L];
        d_181  =  d_180 / 57.29577951308232;
        ul_3[56L]  =  d_181;
        d_182  =  ul_3[64L];
        d_183  =  d_182 / 57.29577951308232;
        ul_3[64L]  =  d_183;
      }  // B20
      else
      {
        // BLOCK 21
      }  // B21
      // BLOCK 22 MERGES [21 20 ]
      b_184  =  i_38 < 2;
      if(b_184)
      {
        // BLOCK 23
        d_185  =  ul_3[l_42];
        d_186  =  ul_3[l_47];
        d_187  =  ul_3[56L];
        d_188  =  ul_3[64L];
        ul_189  =  *((__global ulong *) ul_37);
        ul_190  =  ul_1 + ul_189;
        i_191  =  i_174 + 3;
        l_192  =  (long) i_191;
        l_193  =  l_192 << 3;
        ul_194  =  ul_190 + l_193;
        *((__global double *) ul_194)  =  d_185;
        ul_195  =  *((__global ulong *) ul_37);
        ul_196  =  ul_1 + ul_195;
        i_197  =  i_174 + 4;
        l_198  =  (long) i_197;
        l_199  =  l_198 << 3;
        ul_200  =  ul_196 + l_199;
        *((__global double *) ul_200)  =  d_186;
        ul_201  =  *((__global ulong *) ul_37);
        ul_202  =  ul_1 + ul_201;
        i_203  =  i_174 + 5;
        l_204  =  (long) i_203;
        l_205  =  l_204 << 3;
        ul_206  =  ul_202 + l_205;
        *((__global double *) ul_206)  =  d_187;
        ul_207  =  *((__global ulong *) ul_37);
        ul_208  =  ul_1 + ul_207;
        i_209  =  i_174 + 6;
        l_210  =  (long) i_209;
        l_211  =  l_210 << 3;
        ul_212  =  ul_208 + l_211;
        *((__global double *) ul_212)  =  d_188;
        ul_213  =  *((__global ulong *) ul_37);
        ul_214  =  ul_1 + ul_213;
        ul_215  =  ul_214 + l_177;
        d_216  =  native_sin(d_185);
        d_217  =  native_sin(d_187);
        d_218  =  native_cos(d_185);
        d_219  =  native_cos(d_187);
        d_220  =  d_218 * d_219;
        d_221  =  d_186 - d_188;
        d_222  =  native_cos(d_221);
        d_223  =  d_220 * d_222;
        d_224  =  fma(d_216, d_217, d_223);
        d_225  =  acos(d_224);
        d_226  =  d_225 * 6366707.019493707;
        *((__global double *) ul_215)  =  d_226;
      }  // B23
      else
      {
        // BLOCK 24
        ul_227  =  *((__global ulong *) ul_37);
        ul_228  =  ul_1 + ul_227;
        ul_229  =  ul_228 + l_177;
        *((__global double *) ul_229)  =  0.0;
      }  // B24
    }  // B22
    else
    {
      // BLOCK 25
      ul_230  =  *((__global ulong *) ul_37);
      ul_231  =  ul_1 + ul_230;
      ul_232  =  ul_231 + l_177;
      *((__global double *) ul_232)  =  0.0;
    }  // B25
    // BLOCK 26 MERGES [23 25 24 ]
    b_233  =  i_38 == 0;
    if(b_233)
    {
      // BLOCK 27
      d_234  =  ul_3[72L];
      d_235  =  d_234 / 57.29577951308232;
      ul_3[72L]  =  d_235;
      d_236  =  ul_3[80L];
      d_237  =  d_236 / 57.29577951308232;
      ul_3[80L]  =  d_237;
    }  // B27
    else
    {
      // BLOCK 28
    }  // B28
    // BLOCK 29 MERGES [28 27 ]
    i_238  =  i_38 + 1;
    i_239  =  i_50 + 3;
    i_240  =  i_239 << 2;
    i_241  =  i_240 + i_239;
    i_242  =  i_241 + 7;
    l_243  =  (long) i_242;
    l_244  =  l_243 << 3;
    b_245  =  i_38 < 3;
    if(b_245)
    {
      // BLOCK 30
      d_246  =  ul_3[l_42];
      d_247  =  ul_3[l_47];
      d_248  =  ul_3[72L];
      d_249  =  ul_3[80L];
      ul_250  =  *((__global ulong *) ul_37);
      ul_251  =  ul_1 + ul_250;
      i_252  =  i_241 + 3;
      l_253  =  (long) i_252;
      l_254  =  l_253 << 3;
      ul_255  =  ul_251 + l_254;
      *((__global double *) ul_255)  =  d_246;
      ul_256  =  *((__global ulong *) ul_37);
      ul_257  =  ul_1 + ul_256;
      i_258  =  i_241 + 4;
      l_259  =  (long) i_258;
      l_260  =  l_259 << 3;
      ul_261  =  ul_257 + l_260;
      *((__global double *) ul_261)  =  d_247;
      ul_262  =  *((__global ulong *) ul_37);
      ul_263  =  ul_1 + ul_262;
      i_264  =  i_241 + 5;
      l_265  =  (long) i_264;
      l_266  =  l_265 << 3;
      ul_267  =  ul_263 + l_266;
      *((__global double *) ul_267)  =  d_248;
      ul_268  =  *((__global ulong *) ul_37);
      ul_269  =  ul_1 + ul_268;
      i_270  =  i_241 + 6;
      l_271  =  (long) i_270;
      l_272  =  l_271 << 3;
      ul_273  =  ul_269 + l_272;
      *((__global double *) ul_273)  =  d_249;
      ul_274  =  *((__global ulong *) ul_37);
      ul_275  =  ul_1 + ul_274;
      ul_276  =  ul_275 + l_244;
      d_277  =  native_sin(d_246);
      d_278  =  native_sin(d_248);
      d_279  =  native_cos(d_246);
      d_280  =  native_cos(d_248);
      d_281  =  d_279 * d_280;
      d_282  =  d_247 - d_249;
      d_283  =  native_cos(d_282);
      d_284  =  d_281 * d_283;
      d_285  =  fma(d_277, d_278, d_284);
      d_286  =  acos(d_285);
      d_287  =  d_286 * 6366707.019493707;
      *((__global double *) ul_276)  =  d_287;
    }  // B30
    else
    {
      // BLOCK 31
      ul_288  =  *((__global ulong *) ul_37);
      ul_289  =  ul_1 + ul_288;
      ul_290  =  ul_289 + l_244;
      *((__global double *) ul_290)  =  0.0;
    }  // B31
    // BLOCK 32 MERGES [30 31 ]
    i_291  =  i_238;
    i_38  =  i_291;
  }  // B32
  // BLOCK 33
  return;
}  //  kernel
