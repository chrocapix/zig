// Implementation of:
//
// * "Table-Driven Implementation of the Exponential Function in IEEE Floating-Point Arithmetic"
//   ACM Transactions on Mathematical Software (TOMS), Volume 15, Issue 2 pp. 144-157
//
// and
//
// * "Table-driven implementation of the Expm1 function in IEEE floating-point arithmetic"
//   ACM Transactions on Mathematical Software (TOMS), Volume 18, Issue 2 pp. 211-222
//
// Both by Ping Tak Peter Tang.
//
// Adapted by Christophe Delage to work with 128 bit floats and with base 2 and 10.
//
// Accurate to <= 0.5 ulp in ~98.5% of cases and to <= 0.6 ulp in 99.99% of cases.
// Worst error seen is <= 0.8 ulp.

const std = @import("std");
const math = std.math;

/// Computes e^x
pub fn expq(x: f128) callconv(.c) f128 {
    // __case = 0;
    if (!math.isFinite(x)) {
        if (math.isNan(x)) {
            if (math.isSignalNan(x)) math.raiseInvalid();
            return math.nan(f128);
        }
        return if (math.signbit(x)) 0 else x;
    }
    if (@abs(x) > 12000)
        return if (math.signbit(x)) underflow() else overflow();
    if (@abs(x) < 0x1p-114)
        return 1 + x;

    return proc1(cfg_e, x);
}

/// Computes e^x - 1, with high precision near x = 0.
pub fn expm1q(x: f128) callconv(.c) f128 {
    // __case = 0;
    if (!math.isFinite(x)) {
        if (math.isNan(x)) {
            if (math.isSignalNan(x)) math.raiseInvalid();
            return math.nan(f128);
        }
        return if (math.signbit(x)) -1 else x;
    }
    if (@abs(x) < 0x1p-114) {
        if (@abs(x) <= 0) return 0;
        return math.ldexp(math.ldexp(x, 200) + @abs(x), -200);
    }
    if (x < -79.0187)
        return -1.0 + underflow();
    if (x > 31191.6)
        return overflow();

    const lo = -0.28768207245178092743921900599382744;
    const hi = 0.22314355131420975576629509030983452;
    if (lo < x and x < hi)
        return expm1Proc2(x);

    return proc1(cfg_em1, x);
}

/// Computes 2^x.
pub fn exp2q(x: f128) callconv(.c) f128 {
    // __case = 0;
    if (!math.isFinite(x)) {
        if (math.isNan(x)) {
            if (math.isSignalNan(x)) math.raiseInvalid();
            return math.nan(f128);
        }
        return if (math.signbit(x)) 0 else x;
    }
    if (@abs(x) > 17000)
        return if (math.signbit(x)) underflow() else overflow();
    if (@abs(x) < 0x1p-114 * math.log2e)
        return 1 + x * math.ln2;

    return proc1(cfg_2, x);
}

/// Computes 10^x.
pub fn exp10q(x: f128) callconv(.c) f128 {
    // __case = 0;
    if (!math.isFinite(x)) {
        if (math.isNan(x)) {
            if (math.isSignalNan(x)) math.raiseInvalid();
            return math.nan(f128);
        }
        return if (math.signbit(x)) 0 else x;
    }
    if (@abs(x) > 5000)
        return if (math.signbit(x)) underflow() else overflow();
    if (@abs(x) < 0x1p-114 * math.log10e)
        return 1 + x * math.ln10;

    return proc1(cfg_10, x);
}

fn underflow() f128 {
    const tiny: f128 = 0x1p-16000;
    const p: *const volatile f128 = &tiny;
    return p.* * p.*;
}

fn overflow() f128 {
    const huge: f128 = 0x1p16000;
    const p: *const volatile f128 = &huge;
    return p.* * p.*;
}

const cfg_e: Config = .{
    .inv_log_32 = 0x1.71547652b82fe1777d0ffda0d23ap5,
    .log_32_hi = 0x1.62e42fefa39ef35793c7673p-6,
    .neg_log_32_lo = -0x1.f97b57a079a193394c5b16c5068cp-108,
    .poly = expPoly,
    .finalize = expFinalize,
};

const cfg_em1: Config = .{
    .inv_log_32 = 0x1.71547652b82fe1777d0ffda0d23ap5,
    .log_32_hi = 0x1.62e42fefa39ef35793c7673p-6,
    .neg_log_32_lo = -0x1.f97b57a079a193394c5b16c5068cp-108,
    .poly = expPoly,
    .finalize = expm1Finalize,
};

const cfg_2: Config = .{
    .inv_log_32 = 32.0,
    .log_32_hi = 1.0 / 32.0,
    .neg_log_32_lo = 0.0,
    .poly = exp2Poly,
    .finalize = expFinalize,
};

const cfg_10: Config = .{
    .inv_log_32 = 1.0630169903639559513185022174366049e2,
    .log_32_hi = 0x1.34413509f79fef311f12b358p-7,
    .neg_log_32_lo = -0x1.6f922f04d5a618a87a3e69314bcep-107,
    .poly = exp10Poly,
    .finalize = expFinalize,
};

/// computes (1 + pr) * 2^(j / 32) * 2^m.
fn expFinalize(pr: f128, j: usize, m: i32) f128 {
    const sj_hi = exp2_32[j].hi;
    const sj_lo = exp2_32[j].lo;
    const sj = sj_hi + sj_lo;

    return math.ldexp(sj_hi + (sj_lo + sj * pr), m);
}

/// Approximates e^(r_hi + r_lo) - 1 in [-log(2) / 64, log(2) / 64]
fn expPoly(r_hi: f128, r_lo: f128) f128 {
    // e^r- 1 ~= r + sum(ai * r^i, i = 2..11)
    // theoretical error <= 1.0620373e-35 in [-log(2) / 64, log(2) / 64]
    const a2 = 0.5000000000000000000000000000000905;
    const a3 = 0.1666666666666666666666666666667363;
    const a4 = 4.1666666666666666666666628071963674e-2;
    const a5 = 8.333333333333333333333323239792071e-3;
    const a6 = 1.388888888888888891521136130433685e-3;
    const a7 = 1.9841269841269841311955869505827727e-4;
    const a8 = 2.4801587301524467738593346025740415e-5;
    const a9 = 2.7557319223912701135760993821204538e-6;
    const a10 = 2.7557380444148826192829010473072314e-7;
    const a11 = 2.505216489575622432728510602024213e-8;

    const r = r_hi + r_lo;
    const r2 = r * r;
    const s: f64 = @floatCast(r);

    const a8_9 = a8 + s * a9;
    const a10_11 = a10 + s * a11;
    const a8_11 = a8_9 + s * s * a10_11;
    const a6_7 = a6 + r * a7;
    const a4_5 = a4 + r * a5;
    const a6_11 = a6_7 + r2 * a8_11;
    const a2_3 = a2 + r * a3;
    const a4_11 = a4_5 + r2 * a6_11;
    const a2_11 = a2_3 + r2 * a4_11;
    return r_hi + (r_lo + r2 * a2_11);
}

/// Approximates e^r - 1 in [log(1 - 1 / 4), log(1 + 1 / 4)]
fn expm1Proc2(r: f128) f128 {
    // e^r - 1 ~= r + r^2 / 2 + sum(ai * r^i, i = 3..20)
    // theoretical error <= 8.277959e-39 in [log(1 - 1 / 4), log(1 + 1 / 4)]
    const a3 = 0.1666666666666666666666666666666665;
    const a4 = 4.1666666666666666666666666666664594e-2;
    const a5 = 8.333333333333333333333333333428546e-3;
    const a6 = 1.3888888888888888888888888897107884e-3;
    const a7 = 1.984126984126984126984126809960433e-4;
    const a8 = 2.480158730158730158730145664181578e-5;
    const a9 = 2.7557319223985890652572854614548794e-6;
    const a10 = 2.755731922398589065366458695852114e-7;
    const a11 = 2.5052108385441718699275801593735587e-8;
    const a12 = 2.087675698786809344851457445175285e-9;
    const a13 = 1.6059043836821819849839849769808023e-10;
    const a14 = 1.1470745597746522081495308596851777e-11;
    const a15 = 7.64716373154474168559478905654385e-13;
    const a16 = 4.7794773019257805811249254796805686e-14;
    const a17 = 2.811457318435321051617817752049882e-15;
    const a18 = 1.5619509745996214937038392219702005e-16;
    const a19 = 8.22306549140949601031988140517731e-18;
    const a20 = 3.984099431604597427589421326297223e-19;

    const r2 = r * r;
    const r3 = r2 * r;
    const r4 = r2 * r2;
    const s: f64 = @floatCast(r);

    const a19_20 = a19 + s * a20;
    const a17_18 = a17 + s * a18;
    const a15_16 = a15 + s * a16;
    const a13_14 = a13 + r * a14;
    const a11_12 = a11 + r * a12;
    const a9_10 = a9 + r * a10;
    const a7_8 = a7 + r * a8;
    const a5_6 = a5 + r * a6;
    const a3_4 = a3 + r * a4;
    const a17_20 = a17_18 + s * s * a19_20;
    const a13_16 = a13_14 + r2 * a15_16;
    const a9_12 = a9_10 + r2 * a11_12;
    const a3_8 = a3_4 + r2 * (a5_6 + r2 * a7_8);
    const q = r3 * (a3_8 + r2 * r4 * (a9_12 + r4 * (a13_16 + r4 * a17_20)));

    const u: f128 = @as(f64, @floatCast(r));
    const v = r - u;
    const y = u * u * 0.5;
    const z = v * (r + u) * 0.5;
    return r + (y + (q + z));
}

/// Computes (1 + pr) * 2^(j / 32 + m) - 1
fn expm1Finalize(pr: f128, j: usize, m: i32) f128 {
    const sj_hi = exp2_32[j].hi;
    const sj_lo = exp2_32[j].lo;
    const sj = sj_hi + sj_lo;

    const inv_2m = math.ldexp(@as(f128, 1), -m);

    // Careful computation of (sj * (1 + pr)) * 2^m - 1, depending on the value
    // of m
    if (m >= 113) {
        const frac = sj_hi + ((sj_lo - inv_2m) + sj * pr);
        return math.ldexp(frac, m);
    }
    if (m <= -7) {
        const frac = sj_hi + (sj_lo + sj * pr);
        return math.ldexp(frac, m) - 1;
    }
    const frac = (sj_hi - inv_2m) + (sj_lo * (1 + pr) + sj_hi * pr);
    return math.ldexp(frac, m);
}

/// Approximates 2^(r_hi + r_lo) - 1 in [-1 / 64, 1 / 64]
fn exp2Poly(r_hi: f128, r_lo: f128) f128 {
    // 2^r - 1 ~= sum(ai * r^i, i = 1..11)
    // theoretical error <= 5.3680108e-36 in [-1 / 64, 1 / 64]
    const a1 = 0.6931471805599453094172321214581766;
    const a2 = 0.24022650695910071233355126316357173;
    const a3 = 5.5504108664821579953142263768692476e-2;
    const a4 = 9.618129107628477161979051973829494e-3;
    const a5 = 1.3333558146428443423412192542693712e-3;
    const a6 = 1.5403530393381609999394397511606273e-4;
    const a7 = 1.5252733804059840325781982612177378e-5;
    const a8 = 1.321548679010221918742158270799948e-6;
    const a9 = 1.0178086009207059129274011065740134e-7;
    const a10 = 7.054928860981789930134373497307763e-9;
    const a11 = 4.445549135108216233399948317452361e-10;

    const r = r_hi + r_lo;
    const r2 = r * r;
    const s: f64 = @floatCast(r);

    const a10_11 = a10 + s * a11;
    const a8_9 = a8 + s * a9;
    const a8_11 = a8_9 + s * s * a10_11;
    const a6_7 = a6 + r * a7;
    const a4_5 = a4 + r * a5;
    const a2_3 = a2 + r * a3;
    return r * a1 + r2 * (a2_3 + r2 * (a4_5 + r2 * (a6_7 + r2 * a8_11)));
}

/// Approximates 10^(r_hi + r_lo) - 1 in [-log_10(2) / 64, log_10(2) / 64]
fn exp10Poly(r_hi: f128, r_lo: f128) f128 {
    // 10^r - 1 ~= sum(ai * r^i, i = 1..11)
    // theoretical error <= 5.3680108e-36 in [-log_10(2) / 64, log_10(2) / 64]
    const a1 = 2.302585092994045684017991454684364;
    const a2 = 2.650949055239199005280833194299649;
    const a3 = 2.034678592293476196830991191716403;
    const a4 = 1.1712551489122669631782552292481945;
    const a5 = 0.5393829291955814101996903645582445;
    const a6 = 0.20699584869686809730381157655679514;
    const a7 = 6.808936507443706256966473488963759e-2;
    const a8 = 1.959769462641610655862603447683257e-2;
    const a9 = 5.013928833759360686244531585463725e-3;
    const a10 = 1.1545026002647865053408777040852769e-3;
    const a11 = 2.4166731608758542205843156038924033e-4;

    const r = r_hi + r_lo;
    const r2 = r * r;
    const s: f64 = @floatCast(r);

    const a10_11 = a10 + s * a11;
    const a8_9 = a8 + s * a9;
    const a8_11 = a8_9 + s * s * a10_11;
    const a6_7 = a6 + r * a7;
    const a4_5 = a4 + r * a5;
    const a2_3 = a2 + r * a3;
    return r * a1 + r2 * (a2_3 + r2 * (a4_5 + r2 * (a6_7 + r2 * a8_11)));
}

/// Configuration for computing a^x, a in {e, 2, 10}, or e^x - 1.
const Config = struct {
    /// 32 / log(a)
    inv_log_32: f128,

    // High bits of log(a) / 32
    log_32_hi: f128,
    // Low bits of log(a) / 32, negated
    neg_log_32_lo: f128,

    /// Approximates a^(r_hi + r_lo) - 1.
    poly: fn (r_hi: f128, r_lo: f128) f128,
    /// Compute the final value, a^x or e^x - 1.
    finalize: fn (pr: f128, j: usize, m: i32) f128,
};

/// Compute a^`x` or a^`x` - 1, depending on `cfg`.
fn proc1(comptime cfg: Config, x: f128) f128 {
    // Argument reduction: x = r * 2^(j / 32 + m)
    // with r in [-log_a(2) / 64, log_a(2) / 64].
    //
    // r computed as r_hi + r_lo to simulate higher precision.
    const n = @round(x * cfg.inv_log_32);
    const ni: i32 = @intFromFloat(n);
    // __case = ni;
    const n2 = @mod(ni, 32);
    const n1 = ni - n2;
    const m = @divExact(n1, 32);
    const j: usize = @intCast(n2);

    const r_hi = x - n * cfg.log_32_hi;
    const r_lo = n * cfg.neg_log_32_lo;

    const pr = cfg.poly(r_hi, r_lo);

    return cfg.finalize(pr, j, m);
}

// export var __case: i32 = 0;

/// exp2_32[j] ~= 2^(j / 32)
const exp2_32 = [32]struct { hi: f128, lo: f128 }{
    .{ .hi = 0x1p0, .lo = 0 },
    .{ .hi = 0x1.059b0d31585743ae7c548eb68cap0, .lo = 0x1.05ff94f8d257df7d2ebe128178a8p-110 },
    .{ .hi = 0x1.0b5586cf9890f6298b92b71842ap0, .lo = 0x1.306c8522811679d614545750acd9p-109 },
    .{ .hi = 0x1.11301d0125b50a4ebbf1aed9318p0, .lo = 0x1.9d58b988f562cddcae84e2f0d62p-109 },
    .{ .hi = 0x1.172b83c7d517adcdf7c8c50eb14p0, .lo = 0x1.4f2406aa13fefaeb0d27e4784d7cp-109 },
    .{ .hi = 0x1.1d4873168b9aa7805b8028990fp0, .lo = 0x1.ea62d0881b91859b3c1475c8671cp-110 },
    .{ .hi = 0x1.2387a6e75623866c1fadb1c15cap0, .lo = 0x1.593b0328566902df69e4de240796p-108 },
    .{ .hi = 0x1.29e9df51fdee12c25d15f5a24aap0, .lo = 0x1.de544856046901ff6c05035fb635p-111 },
    .{ .hi = 0x1.306fe0a31b7152de8d5a46305c8p0, .lo = 0x1.7b7b2f09cd0d8a7d40bc6bc5b2e6p-110 },
    .{ .hi = 0x1.371a7373aa9caa7145502f45478p0, .lo = 0x1.87e3e12516bf9c699be432ed3b36p-108 },
    .{ .hi = 0x1.3dea64c12342235b41223e13d76p0, .lo = 0x1.3fba2cb82b8244267c54443f2fp-108 },
    .{ .hi = 0x1.44e086061892d03136f409df018p0, .lo = 0x1.fbd4f3b48709b78591d5cb4f4f32p-108 },
    .{ .hi = 0x1.4bfdad5362a271d4397afec42e2p0, .lo = 0x1.c06c7745c2b38af3f05c9627262ap-113 },
    .{ .hi = 0x1.5342b569d4f81df0a83c49d86a6p0, .lo = 0x1.fa733951f214c02d824a325c9e22p-111 },
    .{ .hi = 0x1.5ab07dd48542958c93015191eb2p0, .lo = 0x1.45d88d7c81280e069fbdb62cbe28p-108 },
    .{ .hi = 0x1.6247eb03a5584b1f0fa06fd2da4p0, .lo = 0x1.5d8e757cfb9913adc57797ced891p-111 },
    .{ .hi = 0x1.6a09e667f3bcc908b2fb1366ea8p0, .lo = 0x1.57d3e3adec17512775099da2f591p-108 },
    .{ .hi = 0x1.71f75e8ec5f73dd2370f2ef0accp0, .lo = 0x1.6cb434b562d9e8a20adda645af12p-108 },
    .{ .hi = 0x1.7a11473eb0186d7d51023f6cdap0, .lo = 0x1.f5ef42b66977960531e821b3497cp-108 },
    .{ .hi = 0x1.82589994cce128acf88afab34ap0, .lo = 0x1.0f6ad65cbbac0f532d39bdfdfdcep-112 },
    .{ .hi = 0x1.8ace5422aa0db5ba7c55a192c9ap0, .lo = 0x1.b3e6ed61f2733304a346d8ed0c01p-108 },
    .{ .hi = 0x1.93737b0cdc5e4f4501c3f2540a2p0, .lo = 0x1.697e257ac0db1f419377f4dd024p-111 },
    .{ .hi = 0x1.9c49182a3f0901c7c46b071f2bep0, .lo = 0x1.6376b7943085c61b242d15de9cc5p-110 },
    .{ .hi = 0x1.a5503b23e255c8b424491caf87ap0, .lo = 0x1.c8050a405381703ef7caff4ecc6p-108 },
    .{ .hi = 0x1.ae89f995ad3ad5e8734d1773204p0, .lo = 0x1.a7fbc3ae675ea440b162d6b8275bp-108 },
    .{ .hi = 0x1.b7f76f2fb5e46eaa7b081ab53c4p0, .lo = 0x1.354c8903c356e4b625aacc2761a4p-108 },
    .{ .hi = 0x1.c199bdd85529c2220cb12a091bap0, .lo = 0x1.99e51125928d998490010e5cceap-110 },
    .{ .hi = 0x1.cb720dcef90691503cbd1e949dap0, .lo = 0x1.761d9559ac0cb6dd3ed599df9ab7p-108 },
    .{ .hi = 0x1.d5818dcfba48725da05aeb66e0cp0, .lo = 0x1.ca9f589f559c0876ff2382fb1a3dp-108 },
    .{ .hi = 0x1.dfc97337b9b5eb968cac39ed29p0, .lo = 0x1.b7225a944efd5bb5524b927250e1p-108 },
    .{ .hi = 0x1.ea4afa2a490d9858f73a18f5db2p0, .lo = 0x1.01f86dea20610ceee13eb7bb0065p-108 },
    .{ .hi = 0x1.f50765b6e4540674f84b762862ap0, .lo = 0x1.aff99000dfc4352ba29b89089653p-108 },
};
