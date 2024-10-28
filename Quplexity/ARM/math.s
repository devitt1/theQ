// © Jacob Liam Gill 2024. All rights reserved. **DO NOT REMOVE THIS LINE.**
// Quplexity MUST be credited in the project you use it in either throughout documentation and/or in the code base. **DO NOT REMOVE THIS LINE.**

// I'm striving to make the code for Quplexity as readable as possible.
// If you would like to contribute or contact me for any other reason please don't hesitate to email me: jacobygill@outlook.com
// Or DM/friend request me on Discord: @mrgill0651

.global _num_pow
.global _esigx01
.global _esigx00
.global _esigx10
//.global _esigx11
.align 3

_num_pow:
    // Expecting nq_L to be in W0
    MOV W1, #1        
    LSL W0, W1, W0    // Shift W1 left by the value in W0 (result = 1 << nq_L)
    RET               

// 0 means first var is pos, 1 means second number is neg (or subtract not add).
// e.g. c * z1 - s * z2
_esigx01:
    // Input: c (D0), s (D1), z1 (D2), z2 (D3)
    FMUL D0, D0, D2     
    FMUL D1, D1, D3     
    FSUB D0, D0, D1     

    RET                 

// e.g. c * z1 + s * z2
_esigx00:
    // Input: c (D0), s (D1), z1 (D2), z2 (D3)
    FMUL D0, D0, D2
    FMUL D1, D1, D3
    FADD D0, D0, D1

    RET

// e.g. -s * z1.y + c * z2.x
_esigx10:
    // Input: -s (D0), c (D1), z1 (D2), z2 (D3)
    FMUL D0, D0, D2
    FMUL D1, D1, D3
    FADD D0, D0, D1

    RET
