//
// Implementation of "Table-driven implementation of the logarithm function in
// IEEE floating-point arithmetic" by Ping-Tak Peter Tang ACM Transactions on
// Mathematical Software (TOMS), Volume 16, Issue 4, pp. 378-400.
//
// Modified by Christophe Delage to work with 128-bit numbers and in base 2 and
// base 10.
//
// Accurate to <= 0.5 ulp in 99.999% of cases. Worst error seen is <= 0.7 ulp.

const std = @import("std");
const math = std.math;

pub fn logq(x: f128) callconv(.c) f128 {
    return log128(config_e, x);
}

pub fn log2q(x: f128) callconv(.c) f128 {
    return log128(config_2, x);
}

pub fn log10q(x: f128) callconv(.c) f128 {
    return log128(config_10, x);
}

pub fn log1pq(x: f128) callconv(.c) f128 {
    // std.debug.print("x = {}\n", .{x});
    if (!math.isFinite(x)) {
        if (math.isNan(x)) {
            if (math.isSignalNan(x)) math.raiseInvalid();
            return math.nan(f128);
        }
        if (!math.signbit(x)) return x;
    }
    if (x <= -1.0) {
        if (x >= -1.0) {
            math.raiseDivByZero();
            return -math.inf(f128);
        }
        math.raiseInvalid();
        return math.nan(f128);
    }
    if (@abs(x) < 0x1p-113) {
        if (@abs(x) <= 0.0) return x;
        return (8 * x - math.floatEpsAt(f128, 0)) * (1.0 / 8.0);
    }

    if (proc2_min1p < x and x < proc2_max1p)
        return proc2(config_e, x);

    const ym = frexp2(if (x < 0x1p115) x + 1 else x);
    const y = ym.significand;
    const m = ym.exponent;
    return proc1(config_e, computeF1p, x, y, m);
}

fn log128(comptime cfg: Config, x: f128) f128 {
    if (!math.isFinite(x)) {
        if (math.isNan(x)) {
            if (math.isSignalNan(x)) math.raiseInvalid();
            return math.nan(f128);
        }
        if (!math.signbit(x)) return x;
    }
    if (x <= 0.0) {
        if (x >= 0.0) {
            math.raiseDivByZero();
            return -math.inf(f128);
        }
        math.raiseInvalid();
        return math.nan(f128);
    }

    if (proc2_min < x and x < proc2_max)
        return proc2(cfg, x - 1.0);

    const ym = frexp2(x);
    const y = ym.significand;
    const m = ym.exponent;
    return proc1(cfg, computeF0, x, y, m);
}

/// exp(-1 / 16) rounded down
const proc2_min = 0.939413062813475786119710824622305;
/// exp(1 / 16) rounded up
const proc2_max = 1.0644944589178594295633905946428897;
/// exp(-1 / 16) - 1 rounded down
const proc2_min1p = -6.0586937186524213880289175377694914e-2;
/// exp(1 / 16) - 1 rounded up
const proc2_max1p = 6.449445891785942956339059464288967e-2;

/// Configuration for computing log_b(x) and log(1 + x)
const Config = struct {
    /// proc1 polynomial:
    /// log_b(1 + f) ~= sum(ai * u^i, i in {1,3,5,7,9,11})
    /// with u = 2 * f / (2 + f) and f in [-1 / 256, 1 / 256].
    a1: f128,
    a3: f128,
    a5: f128,
    a7: f64,
    a9: f64,
    a11: f64,

    /// proc2 polynomial:
    /// log(1 + f) ~= sum(bi * u^i, i in {1,3,5,7,9,11,13,15,17,19})
    /// where 1 + f = (1 + u / 2) / (1 - u / 2), u = 2 * f / (2 + f)
    /// and 1 + f in [e^(-1 / 16), e^(1 / 16)].
    b1_hi: f128,
    b1_lo: f128,
    b3: f128,
    b5: f128,
    b7: f128,
    b9: f128,
    b11: f128,
    b13: f128,
    b15: f64,
    b17: f64,
    b19: f64,

    /// log1p2[j](.hi + .lo) ~= log_b(1 + j / 128),
    /// with the low 11 bits of log1p2.hi set to zero.
    log1p2: [129]struct { hi: f128, lo: f128 },
};

fn computeF0(x: f128, y: f128, m: i32, F: f128) f128 {
    _ = x;
    _ = m;
    return y - F;
}

fn computeF1p(x: f128, y: f128, m: i32, F: f128) f128 {
    if (m <= -2) return y - F;
    const inv2m = math.ldexp(@as(f128, 1.0), -m);
    if (m <= 112) return (inv2m - F) + inv2m * x;
    return (inv2m * x - F) + inv2m;
}

fn proc1(
    comptime cfg: Config,
    comptime computeF: fn (x: f128, y: f128, m: i32, F: f128) f128,
    x: f128,
    y: f128,
    m: i32,
) f128 {
    // let m, F, f such that x = 2^m * (F + f), F = 1 + j / 128,
    // j integer in {0..128}, we compute
    // log_b(x) = m * log_b(2) + log_b(F) + log_b(1 + f / F)
    // where log_b(1 + j / 128) is stored in cfg.log1p2[j]
    // and log_b(1 + f / F) is approximated by
    // P(2 * f / (2 * F + f)) where P(u) is a polynomial
    const F0 = @round(math.ldexp(y, 7));
    const j0: usize = @intFromFloat(F0);
    const j = j0 - 128;
    // std.debug.print("j = {}\n", .{j});
    const F = math.ldexp(F0, -7);
    const f = computeF(x, y, m, F);

    const xm: f128 = @floatFromInt(m);
    const l_hi = xm * cfg.log1p2[128].hi + cfg.log1p2[j].hi;
    const l_lo = xm * cfg.log1p2[128].lo + cfg.log1p2[j].lo;

    const u = (f + f) / (y + F);
    const v = u * u;
    const w: f64 = @floatCast(v);

    const a = cfg.a3 + v * (cfg.a5 + v * (cfg.a7 + w * (cfg.a9 + w * cfg.a11)));
    const q = u * v * a;

    return l_hi + (u * cfg.a1 + (q + l_lo));
}

/// Approximates log(1 + f) for f in [e^(-1 / 16), e^(1 / 16)].
fn proc2(comptime cfg: Config, f: f128) f128 {
    // With u = 2 * f / (2 + f), log(1+f) = log((1 + u / 2) / (1 - u / 2)),
    // which we approximate by a polynomial in u.
    const g = 1 / (2 + f);
    const u = 2 * f * g;
    const v = u * u;
    const w: f64 = @floatCast(v);
    const uv = u * v;

    // we split intermediate values into xxx_hi and xxx_lo to increase precision.
    const q_hi = uv * cfg.b3;
    const q_lo = uv * v * (cfg.b5 + v * (cfg.b7 + v * (cfg.b9 + v * (cfg.b11 +
        v * (cfg.b13 + v * (cfg.b15 + w * (cfg.b17 + w * cfg.b19)))))));

    const f_hi: f128 = @as(f64, @floatCast(f));
    const f_lo: f128 = f - f_hi;

    const u_hi: f128 = @as(f64, @floatCast(u));
    const u_lo: f128 = ((2 * (f - u_hi) - u_hi * f_hi) - u_hi * f_lo) * g;

    // t = u * log_b(e)
    const t_hi = u_hi * cfg.b1_hi;
    const t_lo = u_lo * cfg.b1_hi + u * cfg.b1_lo;

    // y = t + q
    const y_hi = t_hi + q_hi;
    const y_lo = t_lo + (t_hi - y_hi + q_hi) + q_lo;

    return y_hi + y_lo;
}

/// Return (f, k) such that x = f * 2^k and f in [1,2).
/// Asserts that x is finite and positive.
fn frexp2(x: f128) math.Frexp(f128) {
    std.debug.assert(math.isFinite(x));
    std.debug.assert(x > 0.0);

    const bits: u128 = @bitCast(x);
    const uexp: i32 = @intCast(bits >> 112);

    std.debug.assert(uexp >= 0);

    if (uexp == 0) {
        const shift: u7 = @intCast(@clz(bits) - 15);

        const exp = -@as(i32, shift) - 0x3ffe;
        const frac: f128 = @bitCast((bits << shift) | (0x3fff << 112));
        return .{ .significand = frac, .exponent = exp };
    }

    const exp = uexp - 0x3fff;
    const frac: f128 = @bitCast((0x3fff << 112) | ((bits << 16) >> 16));
    return .{ .significand = frac, .exponent = exp };
}

const config_e: Config = .{
    .a1 = 1,
    .a3 = 8.333333333333333333333333333372414e-2,
    .a5 = 1.249999999999999999999963839743372e-2,
    .a7 = 2.2321428571428572515097318595359542e-3,
    .a9 = 4.340277777635300605611118803507141e-4,
    .a11 = 8.877925718782769769445565656611838e-5,

    .b1_hi = 1,
    .b1_lo = 0,
    .b3 = 8.333333333333333333333333333333581e-2,
    .b5 = 1.249999999999999999999999997455655e-2,
    .b7 = 2.2321428571428571428572328745789477e-3,
    .b9 = 4.340277777777777776216500817402857e-4,
    .b11 = 8.87784090909092440759545734146088e-5,
    .b13 = 1.8780048076832339308077858301484127e-5,
    .b15 = 4.069010449774280288178309893970754e-6,
    .b17 = 8.97568550755477160981619052649713e-7,
    .b19 = 2.0165671588771827537210411918018159e-7,

    .log1p2 = .{
        .{ .hi = 0, .lo = 0 },
        .{ .hi = 0x1.fe02a6b106788fc37690392p-8, .lo = -0x1.1ebe96a61269e0832fc62bc4c7ecp-103 },
        .{ .hi = 0x1.fc0a8b0fc03e3cf9eda74d4p-7, .lo = -0x1.0a8552414fc416fc223acca2ebfp-100 },
        .{ .hi = 0x1.7b91b07d5b11aa927f54c72p-6, .lo = -0x1.287fc46561dfab5bc5cceecb4882p-99 },
        .{ .hi = 0x1.f829b0e7833004cf8fc13c8p-6, .lo = -0x1.0dd605151051eb3220ca52e20939p-100 },
        .{ .hi = 0x1.39e87b9febd5fa9015b202ap-5, .lo = 0x1.7229c8d57ae1e11bd3180e8f43acp-98 },
        .{ .hi = 0x1.77458f632dcfc4634f2a1eep-5, .lo = 0x1.2960b1e4dfb80d9544ec6583eb3ap-99 },
        .{ .hi = 0x1.b42dd711971bec28d14c7dap-5, .lo = -0x1.2645ad50c7672fc0eb08d862221dp-102 },
        .{ .hi = 0x1.f0a30c01162a6617cc9716ep-5, .lo = 0x1.d665e2634d1d34cdd6d49659c614p-98 },
        .{ .hi = 0x1.16536eea37ae0e8625c173ep-4, .lo = -0x1.66d0dc92deb7d2ccbd2caa9640cap-99 },
        .{ .hi = 0x1.341d7961bd1d0929983761p-4, .lo = 0x1.344dd408683b296de119e5439841p-98 },
        .{ .hi = 0x1.51b073f06183f69278e686ap-4, .lo = 0x1.7c8ac25e4e3f04de1f086f5cb4b1p-99 },
        .{ .hi = 0x1.6f0d28ae56b4b9be499b9eep-4, .lo = -0x1.cc937e635e7c2135ef0041b8fda2p-97 },
        .{ .hi = 0x1.8c345d6319b20f5acb42a66p-4, .lo = -0x1.254bca8fd9fc1bf283b3b4b8662dp-100 },
        .{ .hi = 0x1.a926d3a4ad563650bd22a9cp-4, .lo = 0x1.d5263cd4fb3f11769cc680ef5589p-99 },
        .{ .hi = 0x1.c5e548f5bc74315d617ef82p-4, .lo = -0x1.3c9d12c4d32a0e47911a6a7dd383p-97 },
        .{ .hi = 0x1.e27076e2af2e5e9ea87ffe2p-4, .lo = -0x1.61eaa246b143bfe80906a822f768p-104 },
        .{ .hi = 0x1.fec9131dbeabaaa2e5199fap-4, .lo = -0x1.9b638802dc3a8afdbfa273d5b1b4p-97 },
        .{ .hi = 0x1.0d77e7cd08e596697717a4p-3, .lo = 0x1.855d1c484cb4fd8d03860ef880c7p-96 },
        .{ .hi = 0x1.1b72ad52f67a029060468e6p-3, .lo = -0x1.94630310e04decd59a27dab257f2p-97 },
        .{ .hi = 0x1.29552f81ff5234c05dc7102p-3, .lo = -0x1.20b2ef60436f8f081d60452c9fc1p-100 },
        .{ .hi = 0x1.371fc201e8f743bcd96c55ep-3, .lo = 0x1.89fcba07cc9b6a9d4a40ca7dd139p-98 },
        .{ .hi = 0x1.44d2b6ccb7d1e67d3d950f8p-3, .lo = 0x1.f878310595f322e9ec1549f016cp-97 },
        .{ .hi = 0x1.526e5e3a1b437a2e401d6e4p-3, .lo = -0x1.a732c9219ce255ddd55e33a13ac2p-98 },
        .{ .hi = 0x1.5ff3070a793d3c873e20a08p-3, .lo = -0x1.bdb88a032aa38caa020c19ce3f53p-96 },
        .{ .hi = 0x1.6d60fe719d21c8d54765c4cp-3, .lo = 0x1.74378df68c9609738dac1eee4e8ep-96 },
        .{ .hi = 0x1.7ab890210d9091be36b2d6ap-3, .lo = 0x1.820191ff8525362042cad5d8c597p-101 },
        .{ .hi = 0x1.87fa06520c910902009017ep-3, .lo = -0x1.b344296aa3ed252052c50575d64dp-98 },
        .{ .hi = 0x1.9525a9cf456b47641307538p-3, .lo = 0x1.712cec4ca0bed3cf7176694713dbp-96 },
        .{ .hi = 0x1.a23bc1fe2b563193711b07ap-3, .lo = 0x1.33184fc7e36cdb1665bcc38099f3p-96 },
        .{ .hi = 0x1.af3c94e80bff2d8ce601938p-3, .lo = -0x1.9851a262591d1968abcd7ad62774p-98 },
        .{ .hi = 0x1.bc286742d8cd629f9ce890ep-3, .lo = 0x1.ea9e1e2c3dca46c83cd6d19e5e9bp-99 },
        .{ .hi = 0x1.c8ff7c79a9a21ac25d81ef2p-3, .lo = 0x1.ff734495c765ea7411adc1b170f1p-96 },
        .{ .hi = 0x1.d5c216b4fbb915b910d65fap-3, .lo = -0x1.95ff1e1c98c2ed4063968ad2333p-96 },
        .{ .hi = 0x1.e27076e2af2e5e9ea87ffe2p-3, .lo = -0x1.61eaa246b143bfe80906a822f768p-103 },
        .{ .hi = 0x1.ef0adcbdc59365218de5438p-3, .lo = -0x1.ef9bd60a5af6797c5b83c5c35a7dp-96 },
        .{ .hi = 0x1.fb9186d5e3e2a8d55466c38p-3, .lo = -0x1.acb4d7db7f4f77549ba6fd09d341p-97 },
        .{ .hi = 0x1.0402594b4d040dae27bd0b6p-2, .lo = -0x1.16a1bbb899f343f105ee37cafa25p-100 },
        .{ .hi = 0x1.0a324e27390e35f73f7a018p-2, .lo = 0x1.c030e28de035afbe0972de36dde4p-96 },
        .{ .hi = 0x1.1058bf9ae4ad5189fa0ab4cp-2, .lo = 0x1.6639f0a672c5cdf89c895a69f18dp-95 },
        .{ .hi = 0x1.1675cababa60e039cc7d572p-2, .lo = -0x1.f23a3ee07cc17dade1f305659e1cp-95 },
        .{ .hi = 0x1.1c898c16999fafbc68e754p-2, .lo = 0x1.c7778677105dd6bd95ebb26c9559p-97 },
        .{ .hi = 0x1.22941fbcf7965a242853da8p-2, .lo = -0x1.4579a168b2a9642c3c4e8c5c0ec5p-95 },
        .{ .hi = 0x1.2895a13de86a35eb49304fcp-2, .lo = 0x1.03962d6a3aacbe5857ab7e2676d6p-98 },
        .{ .hi = 0x1.2e8e2bae11d309c2cc91a86p-2, .lo = -0x1.efc9864244294826ec8756f65987p-95 },
        .{ .hi = 0x1.347dd9a987d54d645674feep-2, .lo = -0x1.9f7846bbd69fd3133c42a01f0d76p-97 },
        .{ .hi = 0x1.3a64c556945e9c72f35cd74p-2, .lo = 0x1.a11beb7a3cee7e029e46e1334dfep-99 },
        .{ .hi = 0x1.404308686a7e3bd0c127df4p-2, .lo = 0x1.8c94c2b20c13ecec8d16d06a791p-95 },
        .{ .hi = 0x1.4618bc21c5ec27d0b7b37b4p-2, .lo = -0x1.871977d446996d9ffa5de361a509p-95 },
        .{ .hi = 0x1.4be5f957778a0db4c9949f6p-2, .lo = 0x1.fb0c8f5ca9a2c339e97829fab5b3p-95 },
        .{ .hi = 0x1.51aad872df82d09c93d60dp-2, .lo = -0x1.5439dc56161950680a3ed41de066p-96 },
        .{ .hi = 0x1.5767717455a6c549ab6ca0ep-2, .lo = -0x1.9f42ff0747cbcce6c0841fd7ceb8p-96 },
        .{ .hi = 0x1.5d1bdbf5809ca508d8e0f72p-2, .lo = -0x1.eea60c7f4b594bd65b44f6b20634p-104 },
        .{ .hi = 0x1.62c82f2b9c7952f6f5f22a6p-2, .lo = 0x1.ca2e7226c55dd257f44b5002c8cdp-102 },
        .{ .hi = 0x1.686c81e9b14aec442be1014p-2, .lo = 0x1.c70d2c96f6f68e19caccb291fd4dp-95 },
        .{ .hi = 0x1.6e08eaa2ba1e38c139318d8p-2, .lo = -0x1.d07a1f9d2a3058cd3f047b933d48p-95 },
        .{ .hi = 0x1.739d7f6bbd0069ce24c53fap-2, .lo = 0x1.a7def2d56ef086e752cbabddde49p-95 },
        .{ .hi = 0x1.792a55fdd47a27c15da47fap-2, .lo = 0x1.ed68159fc32ef1838fd7bd5ea6f5p-96 },
        .{ .hi = 0x1.7eaf83b82afc364b3a5e7b4p-2, .lo = 0x1.4a821b8b065fe5fb8f746a59e683p-95 },
        .{ .hi = 0x1.842d1da1e8b17493b1465e2p-2, .lo = -0x1.cc4cbb87dcf5c7063f8350966c46p-95 },
        .{ .hi = 0x1.89a3386c1425ab5a718811p-2, .lo = 0x1.e0ab67ff8c6503680a2b463ea03fp-97 },
        .{ .hi = 0x1.8f11e873662c77e1769d56ap-2, .lo = -0x1.e5d6686340343f1babef11735f2fp-96 },
        .{ .hi = 0x1.947941c2116faba4cdd147ep-2, .lo = -0x1.df22a2b6b19ed11af82f2c0e6731p-95 },
        .{ .hi = 0x1.99d958117e08acba92eec48p-2, .lo = -0x1.f39022da2471b505dae6ed251384p-96 },
        .{ .hi = 0x1.9f323ecbf984bf2b68d767p-2, .lo = -0x1.7f5bbcfcc16f80265d5c141ea63bp-95 },
        .{ .hi = 0x1.a484090e5bb0a2bfca6b70ep-2, .lo = 0x1.69d25939d6f50c6b635ed1468194p-95 },
        .{ .hi = 0x1.a9cec9a9a08498d484ff53p-2, .lo = -0x1.b579420f701b88ecc4b9dd8b0653p-95 },
        .{ .hi = 0x1.af1293247786b1133844a16p-2, .lo = -0x1.1ebf9e3b2fb68378f9b7aa0ef668p-97 },
        .{ .hi = 0x1.b44f77bcc8f628cbeedaaeap-2, .lo = -0x1.c78710dfac54b2f79fc30eee5577p-96 },
        .{ .hi = 0x1.b9858969310fb598fb14f88p-2, .lo = 0x1.de1cf79039d5a31010282d8264bp-95 },
        .{ .hi = 0x1.beb4d9da71b7bf7861d37acp-2, .lo = -0x1.0f94da4ad3e41effdcc4cd6b6d09p-96 },
        .{ .hi = 0x1.c3dd7a7cdad4d73b3c14b7ap-2, .lo = 0x1.34db05da8c72bec0a04cff3ce95fp-95 },
        .{ .hi = 0x1.c8ff7c79a9a21ac25d81ef2p-2, .lo = 0x1.ff734495c765ea7411adc1b170f1p-95 },
        .{ .hi = 0x1.ce1af0b85f3eb7b7d2bcaaep-2, .lo = -0x1.ecc5b1de7373ff4d090d4dc66787p-95 },
        .{ .hi = 0x1.d32fe7e00ebd561dec8cbecp-2, .lo = -0x1.7240b9d52c802840a469f003b619p-98 },
        .{ .hi = 0x1.d83e7258a2f3e50515ba2ecp-2, .lo = 0x1.288875284a660c123bf5ec01a84p-95 },
        .{ .hi = 0x1.dd46a04c1c4a0bee626a49ep-2, .lo = -0x1.c48f00b057cbd98ca28181536d78p-95 },
        .{ .hi = 0x1.e24881a7c6c261cbd8f4596p-2, .lo = -0x1.6b7393abfc20e89b1e12de69bc58p-95 },
        .{ .hi = 0x1.e744261d68787e37da36f3cp-2, .lo = 0x1.768dda0250c71552f700a89054abp-95 },
        .{ .hi = 0x1.ec399d2468cc0175cee53f4p-2, .lo = -0x1.58d2027bb4681af8a138d7763332p-95 },
        .{ .hi = 0x1.f128f5faf06ecb35c83b114p-2, .lo = -0x1.c61451854db7745eb09c5c270ad2p-95 },
        .{ .hi = 0x1.f6123fa7028ac61456c3cb6p-2, .lo = 0x1.9a0a432a2414cc0b049c0fb73b4ep-95 },
        .{ .hi = 0x1.faf588f78f31ed9afb3e4eap-2, .lo = 0x1.0d7f636a66f1778c26e3db731a5dp-95 },
        .{ .hi = 0x1.ffd2e0857f4985597d0364cp-2, .lo = 0x1.0d7923719bce9f534714c3359d9cp-95 },
        .{ .hi = 0x1.02552a5a5d0fec69c695d7ep-1, .lo = 0x1.cfff8d1b92e2e3130d3aee27a895p-94 },
        .{ .hi = 0x1.04bdf9da926d265fcc1008cp-1, .lo = -0x1.b104ac97b8ddc78779c7d607d871p-94 },
        .{ .hi = 0x1.0723e5c1cdf404e57963892p-1, .lo = -0x1.d22a9352168a2b7bb7bc228c5fd8p-94 },
        .{ .hi = 0x1.0986f4f573520b91fda95p-1, .lo = -0x1.7807f1fbb03b297671d83b0ed087p-94 },
        .{ .hi = 0x1.0be72e4252a82b69897bb34p-1, .lo = -0x1.0019649bc990440ca2c12ee56f6dp-96 },
        .{ .hi = 0x1.0e44985d1cc8bf6eae5de96p-1, .lo = 0x1.2c7c68fdd67b8516d0c2f4242cf2p-94 },
        .{ .hi = 0x1.109f39e2d4c96fde3ec9b0cp-1, .lo = -0x1.0dce836fa706419be06d47b1d77ep-94 },
        .{ .hi = 0x1.12f719593efbc53012319c8p-1, .lo = -0x1.795a5db1362a02bfb7321944cbaep-96 },
        .{ .hi = 0x1.154c3d2f4d5e9a98f33a396p-1, .lo = 0x1.78a02a769d19879a0e78fd77137p-95 },
        .{ .hi = 0x1.179eabbd899a0bfc60e6fap-1, .lo = 0x1.f4b86ac71bec323b0a3e36a42df1p-95 },
        .{ .hi = 0x1.19ee6b467c96ecc5cbdd778p-1, .lo = 0x1.de7aae20f1d40d8586e3a6dfd31ep-97 },
        .{ .hi = 0x1.1c3b81f713c24bc94f8ecep-1, .lo = -0x1.0e881c4059ee4326883918e43f55p-95 },
        .{ .hi = 0x1.1e85f5e7040d03dec59a5f4p-1, .lo = -0x1.c3941a306b616fd12d2a58c0cb6ap-97 },
        .{ .hi = 0x1.20cdcd192ab6d93503d0f76p-1, .lo = -0x1.013db0bb5bf2b76be256c1a2a08bp-95 },
        .{ .hi = 0x1.23130d7bebf4282de368722p-1, .lo = 0x1.9490291c6f37f307d6c31b1816a8p-94 },
        .{ .hi = 0x1.2555bce98f7cb3c043adad2p-1, .lo = -0x1.f863bdee2bd7e4533f945f3dc04dp-94 },
        .{ .hi = 0x1.2795e1289b11aeb783f3dbap-1, .lo = -0x1.2f1c00ff2b60a33daf311082e2f2p-94 },
        .{ .hi = 0x1.29d37fec2b08ac85cd6cba6p-1, .lo = -0x1.9824a2f303e0aa6995ef5598cc5bp-96 },
        .{ .hi = 0x1.2c0e9ed448e8bb97a9c31bap-1, .lo = 0x1.4f312a40a546f842b745196d2d84p-96 },
        .{ .hi = 0x1.2e47436e402684054218692p-1, .lo = 0x1.0d589296794a2964475e2e7179b8p-96 },
        .{ .hi = 0x1.307d7334f10be1fb590a1f6p-1, .lo = -0x1.324c8479fd23fb331fb9915bafefp-94 },
        .{ .hi = 0x1.32b1339121d71320556b67cp-1, .lo = -0x1.ba53e63a917b32e74881ad086531p-94 },
        .{ .hi = 0x1.34e289d9ce1d316eb92d886p-1, .lo = -0x1.8d8a8dace2202c7d4941bef6389bp-96 },
        .{ .hi = 0x1.37117b54747b5c5dd024844p-1, .lo = 0x1.bbf23e2563661b0a05b0d9343044p-94 },
        .{ .hi = 0x1.393e0d3562a19a9c4426036p-1, .lo = 0x1.be30e3d885c5f279c0104c98107ap-94 },
        .{ .hi = 0x1.3b68449fffc22af8edec186p-1, .lo = -0x1.b265f1490c93235372d5c3796b29p-95 },
        .{ .hi = 0x1.3d9026a7156faa404263d0ap-1, .lo = 0x1.a4f3dae02765a552b7c7fe4d82f4p-94 },
        .{ .hi = 0x1.3fb5b84d16f425b4e9d505cp-1, .lo = 0x1.57779b213b3cceb2d2dcad3e5c46p-94 },
        .{ .hi = 0x1.41d8fe84672ae6464bcc2f4p-1, .lo = 0x1.805de54e22437513ab7accba123ep-95 },
        .{ .hi = 0x1.43f9fe2f9ce677a727b9b6p-1, .lo = 0x1.dc32db06eddd30bf9f399846b88fp-94 },
        .{ .hi = 0x1.4618bc21c5ec27d0b7b37b4p-1, .lo = -0x1.871977d446996d9ffa5de361a509p-94 },
        .{ .hi = 0x1.48353d1ea88df73d5e8bb3p-1, .lo = 0x1.fd096183df7cf066ea6b94e4db38p-97 },
        .{ .hi = 0x1.4a4f85db03ebb0227bf47a6p-1, .lo = 0x1.d86d8c5b9bdc64762bd1f9d116fcp-94 },
        .{ .hi = 0x1.4c679afccee39b168ecdd32p-1, .lo = -0x1.e409e624f2776441c70aa69e7243p-95 },
        .{ .hi = 0x1.4e7d811b75bb09cb0985646p-1, .lo = -0x1.ea88f2f3a1435dfb8445d3611b0bp-95 },
        .{ .hi = 0x1.50913cc01686b4bcb3a5b0cp-1, .lo = -0x1.52784b4371849d59d8d13487153ep-94 },
        .{ .hi = 0x1.52a2d265bc5aaee77c8af16p-1, .lo = -0x1.3428ac2097b6f84ad9e32670ae3cp-95 },
        .{ .hi = 0x1.54b2467999497a915428b44p-1, .lo = -0x1.3e86976ba2a9508eb1454a066e92p-96 },
        .{ .hi = 0x1.56bf9d5b3f399411c621736p-1, .lo = 0x1.fe59cdc15631bf5c3509451452d4p-96 },
        .{ .hi = 0x1.58cadb5cd79893092f25d94p-1, .lo = -0x1.ee9eabb0e7f953528eee4e120c41p-94 },
        .{ .hi = 0x1.5ad404c359f2cfb29aaa5fp-1, .lo = 0x1.1cd40845839e04578161ccf77754p-96 },
        .{ .hi = 0x1.5cdb1dc6c17648cf6e3c5d8p-1, .lo = -0x1.ef9d427f6bd4735c53c91bc1e36p-94 },
        .{ .hi = 0x1.5ee02a924167570d6095fd2p-1, .lo = 0x1.081f7a6314f8257273d14984f334p-95 },
        .{ .hi = 0x1.60e32f44788d8ca7c895a0cp-1, .lo = -0x1.626aaf32ba0c72294cd5503d5a8p-94 },
        .{ .hi = 0x1.62e42fefa39ef35793c7674p-1, .lo = -0x1.ff0342542fc32f366359d2749d7dp-94 },
    },
};

const config_2: Config = .{
    .a1 = math.log2e,
    .a3 = 0.12022458674074695061332705675072149,
    .a5 = 1.8033688011112042591998536830294507e-2,
    .a7 = 3.22030143055572204818095463930704e-3,
    .a9 = 6.261697225875019234719395591078697e-4,
    .a11 = 1.280813940786848788109850061222256e-4,

    .b1_hi = 0x1.71547652b82fep0,
    .b1_lo = 0x1.777d0ffda0d23a7d11d6aef551bbp-56,
    .b3 = 0.12022458674074695061332705675016125,
    .b5 = 1.8033688011112042591999058475816515e-2,
    .b7 = 3.2203014305557218914285331735164364e-3,
    .b9 = 6.261697226080570342191671010883619e-4,
    .b11 = 1.280801705334664325639770412440281e-4,
    .b13 = 2.7093882228102330360125035037716968e-5,
    .b15 = 5.870341197214724685339102193694838e-6,
    .b17 = 1.294917697032820750200161813672143e-6,
    .b19 = 2.909291439731657940692470637735429e-7,

    .log1p2 = .{
        .{ .hi = 0, .lo = 0 },
        .{ .hi = 0x1.6fe50b6ef08517f8e37b002p-7, .lo = -0x1.0d61777c664136e1b34b7c350515p-100 },
        .{ .hi = 0x1.6e79685c2d2298a6e27e212p-6, .lo = -0x1.fbd41ae7d5a2434912ad3fe21cfbp-100 },
        .{ .hi = 0x1.11cd1d5133412ed814504fap-5, .lo = 0x1.268ea5e6d5ffba8a5f4c2e653545p-98 },
        .{ .hi = 0x1.6bad3758efd87313606f096p-5, .lo = 0x1.6f8212418d9bc96b8da0e67018f5p-98 },
        .{ .hi = 0x1.c4dfab90aab5ef4f8f869e6p-5, .lo = 0x1.dc4142bd2b182fdaac375356fc29p-100 },
        .{ .hi = 0x1.0eb389fa29f9ab3cf74babap-4, .lo = -0x1.9b7a3e651babc30e088b772ba106p-98 },
        .{ .hi = 0x1.3aa2fdd27f1c2d804d1121cp-4, .lo = -0x1.6b3b12b256e5e838744874088ae9p-97 },
        .{ .hi = 0x1.663f6fac913167ccc538262p-4, .lo = -0x1.77514a786b5f060a9c5b00cf6db7p-97 },
        .{ .hi = 0x1.918a16e46335aae7232494ep-4, .lo = -0x1.3171737dca6ec647841a3ca130e1p-98 },
        .{ .hi = 0x1.bc84240adabba63b2c5a6e6p-4, .lo = -0x1.cd0a8f0d37c75efde388fbe6acacp-97 },
        .{ .hi = 0x1.e72ec117fa5b21cbdb5d9dcp-4, .lo = 0x1.4f902752a1dc5b384b68c4f1e669p-99 },
        .{ .hi = 0x1.08c588cda79e39627bc6fdp-3, .lo = 0x1.4aa5157a9e77fd40e0401d0502b1p-96 },
        .{ .hi = 0x1.1dcd197552b7b5ea4543078p-3, .lo = 0x1.cb0abad43eea455a0d24eaa188e2p-98 },
        .{ .hi = 0x1.32ae9e278ae1a1f51f2c076p-3, .lo = -0x1.62d166c9ab03b4ece445b9135432p-97 },
        .{ .hi = 0x1.476a9f983f74d3138e94164p-3, .lo = 0x1.fb609f756d9e3e274454008c86b5p-98 },
        .{ .hi = 0x1.5c01a39fbd6879fa00b120ap-3, .lo = 0x1.a2eb74493cf9a8e8966c101ef964p-101 },
        .{ .hi = 0x1.70742d4ef027f29c01cfad8p-3, .lo = -0x1.03092bad94765b98f2834bbe991cp-96 },
        .{ .hi = 0x1.84c2bd02f03b2fdd2248ee8p-3, .lo = -0x1.38ad51e10526896a04804659b614p-96 },
        .{ .hi = 0x1.98edd077e70df02face8caap-3, .lo = -0x1.d1a901a0697ececc5ac3802a0c6ep-96 },
        .{ .hi = 0x1.acf5e2db4ec93efe11ecbcp-3, .lo = 0x1.83634b52082beea143f1178aaf62p-99 },
        .{ .hi = 0x1.c0db6cdd94dee40e26d989ap-3, .lo = -0x1.8085e509f6df5192d6e95662ccf9p-98 },
        .{ .hi = 0x1.d49ee4c32596fc8f4b56502p-3, .lo = 0x1.d02a20dfc15a5a557fdc6ad4266cp-98 },
        .{ .hi = 0x1.e840be74e6a4cc7c9f3d51ep-3, .lo = -0x1.03c7c0d8554c2c17f944fa1f93cfp-99 },
        .{ .hi = 0x1.fbc16b902680a23a8d998a8p-3, .lo = -0x1.ec432836cc2a4494ad5be52a9827p-97 },
        .{ .hi = 0x1.0790adbb030096f031a699ep-2, .lo = -0x1.5748a4ccf5e3dadf04bd0ff35f88p-95 },
        .{ .hi = 0x1.11307dad30b75cb09705a7ap-2, .lo = -0x1.4b62c07cc46327805ec0076cdce9p-95 },
        .{ .hi = 0x1.1ac05b291f070528c7386ep-2, .lo = -0x1.cd797e78938cb634b4c042530069p-96 },
        .{ .hi = 0x1.24407ab0e07398245b94ba4p-2, .lo = 0x1.300efeef41e86a7632602708133bp-96 },
        .{ .hi = 0x1.2db10fc4d9aaf6f137a3d8cp-2, .lo = 0x1.bcf371dc716701717c89c120ab1bp-96 },
        .{ .hi = 0x1.37124cea4cdecd991336c96p-2, .lo = 0x1.fb9186e41376f4612d4b0b9a507dp-100 },
        .{ .hi = 0x1.406463b1b044975b2f34426p-2, .lo = -0x1.b995b4218dff0ebab8d36942313bp-95 },
        .{ .hi = 0x1.49a784bcd1b8afe492bf7p-2, .lo = -0x1.64a049664d2754039098a562dc93p-95 },
        .{ .hi = 0x1.52dbdfc4c96b37dcf60e62p-2, .lo = -0x1.b2569a211021e1b525598cd22ccbp-97 },
        .{ .hi = 0x1.5c01a39fbd6879fa00b120ap-2, .lo = 0x1.a2eb74493cf9a8e8966c101ef964p-100 },
        .{ .hi = 0x1.6518fe4677ba6e52278edc8p-2, .lo = 0x1.3d310e968bf16ae633ed1b5fa71bp-95 },
        .{ .hi = 0x1.6e221cd9d0cde578d520b44p-2, .lo = 0x1.ee0788864fc48281942fb6c0d1ecp-95 },
        .{ .hi = 0x1.771d2ba7efb3be46fecd512p-2, .lo = 0x1.253fe387923094a7846dddac9a3p-97 },
        .{ .hi = 0x1.800a563161c5432aeb609f4p-2, .lo = 0x1.bef5911b0bd8dfc924d470969c2ep-95 },
        .{ .hi = 0x1.88e9c72e0b225a4b664a4c8p-2, .lo = 0x1.b5375621f63b6e3f997e33b7cec5p-95 },
        .{ .hi = 0x1.91bba891f1708b4b2b5056cp-2, .lo = -0x1.e58ea9e7a24574145c34e7f4c844p-96 },
        .{ .hi = 0x1.9a802391e232f34bb6d0e44p-2, .lo = -0x1.acf95dc3c616f85adff31270d155p-96 },
        .{ .hi = 0x1.a33760a7f60509d7c40d798p-2, .lo = -0x1.84e93808cffe2b0c6f5ecb2d7001p-96 },
        .{ .hi = 0x1.abe18797f1f48e1a4725558p-2, .lo = 0x1.8b7ae3c371e9557e7702b6f5d5dcp-95 },
        .{ .hi = 0x1.b47ebf73882a0a4146ef8fep-2, .lo = -0x1.d761c548eb57fa11bd8e6803cc37p-96 },
        .{ .hi = 0x1.bd0f2e9e79030ab442ce32p-2, .lo = 0x1.16c28a69117eb6d983e10847dac2p-98 },
        .{ .hi = 0x1.c592fad295b567e7ee54aeep-2, .lo = 0x1.e6618a7815661988e883d7a3cdf8p-95 },
        .{ .hi = 0x1.ce0a4923a587cc95d0a2ee8p-2, .lo = -0x1.7e9b7d51fbd6aac640f12dc0b674p-96 },
        .{ .hi = 0x1.d6753e032ea0efe3ebe199p-2, .lo = 0x1.5554d6bf3e730bb7410e895b8a58p-96 },
        .{ .hi = 0x1.ded3fd442364c4ebb196116p-2, .lo = -0x1.df9689b34ee848b0b7818878492cp-99 },
        .{ .hi = 0x1.e726aa1e754d20c519e12f4p-2, .lo = 0x1.c44ccd5e8a375ea75157976bd20ap-96 },
        .{ .hi = 0x1.ef6d67328e2207d1e01a83ap-2, .lo = -0x1.f7db7c031991256093b9dc1c1436p-95 },
        .{ .hi = 0x1.f7a8568cb06cece19318004p-2, .lo = 0x1.43d6c8d5af992540238215bb3ea6p-96 },
        .{ .hi = 0x1.ffd799a83ff9ab9cc7f343p-2, .lo = -0x1.d1e774e84b549ee6e6f32ef31193p-96 },
        .{ .hi = 0x1.03fda8b97997f3394346406p-1, .lo = -0x1.30a5c1a6aa9135b95ba2085b6b35p-94 },
        .{ .hi = 0x1.0809cf27f703d525b3c1d16p-1, .lo = -0x1.f0570f4799a78444f0d0032363d6p-95 },
        .{ .hi = 0x1.0c10500d63aa6588257529cp-1, .lo = -0x1.3b443d55f035e588a5078b8cf48fp-94 },
        .{ .hi = 0x1.10113b153c8ea7b1cddae7p-1, .lo = -0x1.4c6a14b12ca495197c5993ea319dp-95 },
        .{ .hi = 0x1.140c9faa1e5439e15a52a32p-1, .lo = -0x1.3f6b4f76b51289e27c4d51d332d8p-94 },
        .{ .hi = 0x1.18028cf72976a4eb8e97d14p-1, .lo = 0x1.4cd618c18461c458f9708479e9bep-95 },
        .{ .hi = 0x1.1bf311e95d00de3b513a9dcp-1, .lo = 0x1.9b0dc038fab1aaacb83be299f995p-94 },
        .{ .hi = 0x1.1fde3d30e812642415d4738p-1, .lo = 0x1.13c458dd53d12c99743f3c4617c3p-95 },
        .{ .hi = 0x1.23c41d42727c8080ecc61aap-1, .lo = -0x1.dfb113740031e528bbef9ead829cp-95 },
        .{ .hi = 0x1.27a4c0585cbf805784ee0e4p-1, .lo = -0x1.4583145703d2c5dcbd9ccc4a85fbp-95 },
        .{ .hi = 0x1.2b803473f7ad0f3f4016242p-1, .lo = -0x1.7e5d148bb6c30657176993efe107p-94 },
        .{ .hi = 0x1.2f56875eb3f2614278cd16ap-1, .lo = -0x1.b3b9d166fe5d876b1e73c7de470ep-95 },
        .{ .hi = 0x1.3327c6ab49ca6c86b9205fap-1, .lo = 0x1.0880849cf376e0305fc022ceec4p-95 },
        .{ .hi = 0x1.36f3ffb6d9162404772a152p-1, .lo = -0x1.993193dd58663d90eed123bda5eap-96 },
        .{ .hi = 0x1.3abb3faa02166cccab240eap-1, .lo = -0x1.f60d2bdc639cad5bff474c1aa215p-94 },
        .{ .hi = 0x1.3e7d9379f70166ae2a7ada6p-1, .lo = -0x1.5a725949f9fa10048846d5a737c1p-94 },
        .{ .hi = 0x1.423b07e986aa9670761d14ap-1, .lo = 0x1.589b0c986216b63fa1707a772089p-94 },
        .{ .hi = 0x1.45f3a98a20738a4d7ffe026p-1, .lo = 0x1.cbaaa0e3badf44f8420d741bce35p-95 },
        .{ .hi = 0x1.49a784bcd1b8afe492bf7p-1, .lo = -0x1.64a049664d2754039098a562dc93p-94 },
        .{ .hi = 0x1.4d56a5b33cec44a6deff998p-1, .lo = 0x1.cfd68f1bef047af0101b693ac9d2p-95 },
        .{ .hi = 0x1.510118708a8f8dde949378cp-1, .lo = -0x1.bb2dc3aeac5cdb1967f6b532ce4ep-94 },
        .{ .hi = 0x1.54a6e8ca5438db1b0ca63aap-1, .lo = 0x1.768743d98fa7bc6d8c07d9b476f8p-94 },
        .{ .hi = 0x1.5848226989d33c38d8bd28ep-1, .lo = -0x1.2d473ddac42ee7186af4f409cebep-94 },
        .{ .hi = 0x1.5be4d0cb51434aaeb3f0122p-1, .lo = 0x1.12ce7e40053a5cfc072e22bbeab2p-96 },
        .{ .hi = 0x1.5f7cff41e09aeb8cb1ac05cp-1, .lo = 0x1.9ac134ed7c8d8213475a3f8225ffp-94 },
        .{ .hi = 0x1.6310b8f553048406a5a171ep-1, .lo = -0x1.bff3336aeddf91b69ed59b65b4e1p-97 },
        .{ .hi = 0x1.66a008e4788cbcd2edb439p-1, .lo = 0x1.ca60d447873d20f4a4019db1e5d5p-94 },
        .{ .hi = 0x1.6a2af9e5a0f0a08099572f2p-1, .lo = 0x1.7ccd0a8f6177a5b3a9825b18f91p-98 },
        .{ .hi = 0x1.6db196a761949d97df07e36p-1, .lo = -0x1.2bb3cf2d0f250706df598caea05ap-94 },
        .{ .hi = 0x1.7133e9b156c7be5167fbdc8p-1, .lo = 0x1.3726bce2d8f96429416e31f860b1p-97 },
        .{ .hi = 0x1.74b1fd64e0753c6e5783fd2p-1, .lo = -0x1.5e6dc585103c8106450020b13797p-94 },
        .{ .hi = 0x1.782bdbfdda6577bc87e125ep-1, .lo = 0x1.5548be9b13dc830b4c0db4edb473p-94 },
        .{ .hi = 0x1.7ba18f93502e409eab77f22p-1, .lo = -0x1.8e04bf524805be54f509252b5064p-95 },
        .{ .hi = 0x1.7f1322182cf15d12ecd77fep-1, .lo = 0x1.afaed3f53d2caef390a0269b3d52p-95 },
        .{ .hi = 0x1.82809d5be7072dbdc0426c4p-1, .lo = -0x1.ec5cf68c91244c1518f2ef3e8c4fp-96 },
        .{ .hi = 0x1.85ea0b0b27b261086fce864p-1, .lo = 0x1.43eb2f26cfe225c69b789c0f5d25p-94 },
        .{ .hi = 0x1.894f74b06ef8b406ea2c7dap-1, .lo = -0x1.caa76f41d000c510a34d0f88c88dp-94 },
        .{ .hi = 0x1.8cb0e3b4b3bbdb3688a85fcp-1, .lo = -0x1.cc886903f80ca8b19889a5864e8cp-94 },
        .{ .hi = 0x1.900e6160002ccfe43f50848p-1, .lo = -0x1.73ebbc0d5bdb9c4fe022c16a2ecep-96 },
        .{ .hi = 0x1.9367f6da0ab2e9cc865b3dep-1, .lo = -0x1.e4895507862ebfccaf89c26807e2p-94 },
        .{ .hi = 0x1.96bdad2acb5f5efec491532p-1, .lo = -0x1.7241fa312681a6e20d4a0b338b8cp-94 },
        .{ .hi = 0x1.9a0f8d3b0e04fde95734abep-1, .lo = -0x1.a067c67389f4b728e5e6eec32a4ep-94 },
        .{ .hi = 0x1.9d5d9fd5010b36665592074p-1, .lo = 0x1.04f96a11ce31a952005c59f2aef6p-94 },
        .{ .hi = 0x1.a0a7eda4c112ce6312ebb82p-1, .lo = -0x1.8569c9f6eab582c6fdaf350a7a7fp-96 },
        .{ .hi = 0x1.a3ee7f38e181ed0798d1aa2p-1, .lo = 0x1.69456132771589eb06ae62710a0dp-97 },
        .{ .hi = 0x1.a7315d02f20c7bd560a3feep-1, .lo = 0x1.129a29fe9a98fe5dffe5f86e7311p-98 },
        .{ .hi = 0x1.aa708f58014d37cde37c86cp-1, .lo = -0x1.bbe5799290484ee9ad151095340bp-94 },
        .{ .hi = 0x1.adac1e711c832d1562d61bp-1, .lo = -0x1.181015ada3d118153e20dfd4db87p-94 },
        .{ .hi = 0x1.b0e4126bcc86bd7a6ed4e1cp-1, .lo = -0x1.ed94620916c98f8a6d53de3e62d8p-94 },
        .{ .hi = 0x1.b418734a9008bd978b98f7ep-1, .lo = -0x1.2073a650c7a4d0ba364cbd327bd9p-97 },
        .{ .hi = 0x1.b74948f5532da4b4b714336p-1, .lo = 0x1.c637671f05d84c6e2eadff072009p-96 },
        .{ .hi = 0x1.ba769b39e49640ef87ede14p-1, .lo = 0x1.5d1941951c9a376a454e30f028d9p-94 },
        .{ .hi = 0x1.bda071cc67e6db516de0814p-1, .lo = -0x1.325954cfe648ebb98ad18a2ddc38p-94 },
        .{ .hi = 0x1.c0c6d447c5dd362d9a9a55cp-1, .lo = 0x1.d17b370ba83c0155dfdf1fd11697p-95 },
        .{ .hi = 0x1.c3e9ca2e1a05533698b4e4ap-1, .lo = -0x1.213f3f83c76877dcdca4f0a7c286p-95 },
        .{ .hi = 0x1.c7095ae91e1c760bc9b188cp-1, .lo = 0x1.1322631fb315aaf4da97307d1077p-95 },
        .{ .hi = 0x1.ca258dca9331635fee390cp-1, .lo = 0x1.560f40c2c0c5c890acd0f9d8c13fp-94 },
        .{ .hi = 0x1.cd3e6a0ca8906c243749114p-1, .lo = 0x1.895dd7498f46d5338710a594c559p-94 },
        .{ .hi = 0x1.d053f6d2608967318975dcp-1, .lo = 0x1.cf52c6c122a94fa7204a195eb0bp-94 },
        .{ .hi = 0x1.d3663b27f31d5297837adb4p-1, .lo = 0x1.5b1c9574a6344eadef665d73b5ebp-94 },
        .{ .hi = 0x1.d6753e032ea0efe3ebe199p-1, .lo = 0x1.5554d6bf3e730bb7410e895b8a58p-95 },
        .{ .hi = 0x1.d9810643d6614c3c406eb46p-1, .lo = 0x1.105d328adc61c09915e038a135bep-95 },
        .{ .hi = 0x1.dc899ab3ff56c5e673abad4p-1, .lo = 0x1.0c6319cfd3de88693e0324bc5e38p-95 },
        .{ .hi = 0x1.df8f02086af2c4bef483c68p-1, .lo = 0x1.57a2af7075cfa8a2386c6133efa2p-94 },
        .{ .hi = 0x1.e29142e0e01401fbaaa67e4p-1, .lo = -0x1.0eb2a0911dc18a2c9aa3711ef022p-95 },
        .{ .hi = 0x1.e59063c8822ce561911a9bap-1, .lo = 0x1.8cd86f40adb7d8620e56210c92d1p-94 },
        .{ .hi = 0x1.e88c6b3626a72aa21a3c7fp-1, .lo = 0x1.f78e28a80d83e3a4d8e210e57178p-97 },
        .{ .hi = 0x1.eb855f8ca88fb0d4b5c673cp-1, .lo = -0x1.3dc497fc61ad2ce509feb74925dfp-95 },
        .{ .hi = 0x1.ee7b471b3a9507d6dc1f27ep-1, .lo = 0x1.e90f91e68c4501faa71b07ebb524p-94 },
        .{ .hi = 0x1.f16e281db76303b21928c22p-1, .lo = -0x1.412854f795db443461c5e7233e54p-94 },
        .{ .hi = 0x1.f45e08bcf06554e4d5be4f8p-1, .lo = -0x1.f7c0bf059f54635c2e76cded1318p-94 },
        .{ .hi = 0x1.f74aef0efafadd7a1b65f64p-1, .lo = -0x1.cbde0f4c41324535a987d26c23f2p-95 },
        .{ .hi = 0x1.fa34e1177c23362928b9ed8p-1, .lo = -0x1.630adb856a535ad317fb43c83713p-94 },
        .{ .hi = 0x1.fd1be4c7f2af942b221ce0ep-1, .lo = -0x1.df97628deac2911a334146bb3d07p-94 },
        .{ .hi = 0x1p0, .lo = 0 },
    },
};

const config_10: Config = .{
    .a1 = math.log10e,
    .a3 = 3.619120682527098563759407657655348e-2,
    .a5 = 5.428681023790647845638954444458386e-3,
    .a7 = 9.694073256769014481942040422515466e-4,
    .a9 = 1.884958688754320118955531917460363e-4,
    .a11 = 3.8556341504143175800053507804546873e-5,

    .b1_hi = 0x1.bcb7b1526e50ep-2,
    .b1_lo = 0x1.95355baaafad33dc323ee3460246p-57,
    .b3 = 3.619120682527098563759407657638483e-2,
    .b5 = 5.428681023790647845639111475407614e-3,
    .b7 = 9.694073256769014010070232880860484e-4,
    .b9 = 1.8849586888161971679466375197398104e-4,
    .b11 = 3.855597318033137224139238937796866e-5,
    .b13 = 8.156071249646061672592451832809541e-6,
    .b15 = 1.7671487851436387503930903139882346e-6,
    .b17 = 3.898090687230025454479255130305971e-7,
    .b19 = 8.757839894876785986064901881670424e-8,

    .log1p2 = .{
        .{ .hi = 0, .lo = 0 },
        .{ .hi = 0x1.bafd47221ed2665c1ba949p-9, .lo = -0x1.eb6f20a90ad48515635f3b8a1d22p-104 },
        .{ .hi = 0x1.b9476a4fcd10ed89b5a4172p-8, .lo = -0x1.e9d586d6805b5b0784639c343bb5p-101 },
        .{ .hi = 0x1.49b0851443683ce1bf0b25ep-7, .lo = 0x1.7a79983539059615ec7f118f1cf8p-101 },
        .{ .hi = 0x1.b5e908eb137900f974ff1b4p-7, .lo = -0x1.6bbc5e42b470cbc0ccc513d93932p-102 },
        .{ .hi = 0x1.10a83a8446c77a1180aaf6p-6, .lo = -0x1.09ecab3c0cf83110cf07e44c289fp-101 },
        .{ .hi = 0x1.45f4f5acb8be07769e25e96p-6, .lo = -0x1.0a58d59387d112136c57591afc6cp-99 },
        .{ .hi = 0x1.7adc3df3b1ff81b980714c6p-6, .lo = -0x1.a572fb60daf72fbfab402e472463p-100 },
        .{ .hi = 0x1.af5f92b00e60fa6de0a6da8p-6, .lo = -0x1.8335cd731aa841e5f6959299c961p-100 },
        .{ .hi = 0x1.e3806acbd058f0d79f59da2p-6, .lo = 0x1.9ab1b88f11b4d4ea94e47e000ae6p-102 },
        .{ .hi = 0x1.0ba01a81700002be3a8a48ap-5, .lo = -0x1.792eb9b5df590083e09dc0929a7p-101 },
        .{ .hi = 0x1.25502c0fc314b801dad7deep-5, .lo = -0x1.7341c9bdb87272973815e80e69a7p-98 },
        .{ .hi = 0x1.3ed1199a5e425037527d748p-5, .lo = -0x1.c78118e73744326912f10c802895p-100 },
        .{ .hi = 0x1.58238eeb353da7bf5153dfap-5, .lo = 0x1.34ce089c65c3b0105027dcb2ac09p-98 },
        .{ .hi = 0x1.71483427d2a98ce1e11006p-5, .lo = 0x1.c4ea21ef8f4a97aa1e6e70df2d95p-101 },
        .{ .hi = 0x1.8a3fadeb847f393aed3e7b6p-5, .lo = 0x1.098634c8ec2a09971ad607875b6ep-99 },
        .{ .hi = 0x1.a30a9d609efe9c281982d7ep-5, .lo = -0x1.0a32db4c5564f07600c430a0057dp-102 },
        .{ .hi = 0x1.bba9a058dfd841a9796c344p-5, .lo = 0x1.1311440fb8b0d0393066944a085cp-98 },
        .{ .hi = 0x1.d41d5164facb3a0188eb21p-5, .lo = 0x1.0128f682aeb0fa34321a7e97187fp-101 },
        .{ .hi = 0x1.ec6647eb5880847d0188c2cp-5, .lo = 0x1.a9413fcbf87e394c382ab2f6f255p-98 },
        .{ .hi = 0x1.02428c1f08015ea6bc2bc8cp-4, .lo = 0x1.5daed7d59e78eef0247b7ffb956ep-99 },
        .{ .hi = 0x1.0e3d29d81165e62559618f2p-4, .lo = 0x1.9e819128f90a0d3db1dc752ed073p-99 },
        .{ .hi = 0x1.1a23445501815c0cde7a7fp-4, .lo = 0x1.625b060ba6895e77de5083bef908p-98 },
        .{ .hi = 0x1.25f5215eb5949df2a5fb46cp-4, .lo = -0x1.d21e12a99a1833c1b22c98a7c80cp-98 },
        .{ .hi = 0x1.31b3055c4711801b420b9b2p-4, .lo = 0x1.76e9002c845c3b127df25ba8db5bp-103 },
        .{ .hi = 0x1.3d5d335c53178caf84eb228p-4, .lo = 0x1.c747af7739ffd03f91272c27c6a8p-97 },
        .{ .hi = 0x1.48f3ed1df48fb5e08483b68p-4, .lo = -0x1.7f13d7d43d84e36e36265eaea7bcp-101 },
        .{ .hi = 0x1.5477731973e848790c13ee2p-4, .lo = 0x1.0a2a83be3e0b73310670983c73e6p-98 },
        .{ .hi = 0x1.5fe80488af4fca9254c0c64p-4, .lo = -0x1.75244d6c23fa3d7fc5d4ee08bf0dp-98 },
        .{ .hi = 0x1.6b45df6f3e2c9590e0d54c8p-4, .lo = -0x1.3892b6b806179de2dc9f17625c33p-97 },
        .{ .hi = 0x1.769140a2526fc94ecf23d48p-4, .lo = -0x1.e0e7be0fe858091f749bec4fd25ap-97 },
        .{ .hi = 0x1.81ca63d05a449827184d3fep-4, .lo = -0x1.41cca11e23c4b51da2e2aba37379p-97 },
        .{ .hi = 0x1.8cf183886480c9b28b1f97ep-4, .lo = -0x1.13ada343a03c927e5baecee4f8dap-100 },
        .{ .hi = 0x1.9806d9414a2097207328f76p-4, .lo = 0x1.b7b373b2e7ab503dc9f3cd556b4ep-98 },
        .{ .hi = 0x1.a30a9d609efe9c281982d7ep-4, .lo = -0x1.0a32db4c5564f07600c430a0057dp-101 },
        .{ .hi = 0x1.adfd07416be06fd76ea69eep-4, .lo = 0x1.da1b15eecd2991bc415de0639fdbp-98 },
        .{ .hi = 0x1.b8de4d3ab3d97f5dc97fa3ep-4, .lo = 0x1.652bb0f11e7a6d85ba9757d96319p-97 },
        .{ .hi = 0x1.c3aea4a5c6efe9d1b9bf7b4p-4, .lo = 0x1.a458d14aeca9cea45c47fb092e36p-98 },
        .{ .hi = 0x1.ce6e41e463da4f487cfe37cp-4, .lo = -0x1.a180bcb475a3bc1798c9936d0ed7p-97 },
        .{ .hi = 0x1.d91d5866aa99b8c5ecd8544p-4, .lo = 0x1.c860aca4e1162e09ec983aabbdfap-98 },
        .{ .hi = 0x1.e3bc1ab0e19fe3d562a53fp-4, .lo = 0x1.d6361d34de402146e5fb425dbe44p-98 },
        .{ .hi = 0x1.ee4aba610f2047109cc0208p-4, .lo = 0x1.265a92dc7d61d407d46266f8823dp-98 },
        .{ .hi = 0x1.f8c968346819084e03494e8p-4, .lo = -0x1.4b71b85b5d726a32292230bf611dp-99 },
        .{ .hi = 0x1.019c2a064b486717a766838p-3, .lo = 0x1.fd8a0d264f66496899f1a7e9a342p-97 },
        .{ .hi = 0x1.06cbd67a6c3b65458c50fd8p-3, .lo = 0x1.840d79fca5a0c693c952fc0c42bap-103 },
        .{ .hi = 0x1.0bf3d0937c41c3c2f40d06ep-3, .lo = -0x1.24e64860fd68286030ad078843ffp-97 },
        .{ .hi = 0x1.11142f0811356e473b0e4f8p-3, .lo = -0x1.2936162066cf80796938798e0bcfp-96 },
        .{ .hi = 0x1.162d082ac9d0f8e71a2f29p-3, .lo = 0x1.d24d143d5287f3b189160249ade9p-96 },
        .{ .hi = 0x1.1b3e71ec94f7abbbb3324ap-3, .lo = 0x1.ff44989b978452e339c2f0e43dbep-96 },
        .{ .hi = 0x1.204881dee8777552c136a76p-3, .lo = -0x1.267ee5ddb2e8b583793d15630abcp-102 },
        .{ .hi = 0x1.254b4d35e7d3c1d7958ffeep-3, .lo = 0x1.caafb82b7986316a6b6463b79dfcp-97 },
        .{ .hi = 0x1.2a46e8ca7ba29955cdd7838p-3, .lo = 0x1.d028d2ee87e430ccbc0c3839071dp-96 },
        .{ .hi = 0x1.2f3b691c5a000be34bf081ep-3, .lo = 0x1.d538b4bd294d462c8484049e60b9p-97 },
        .{ .hi = 0x1.3428e2540096d3b633e0462p-3, .lo = -0x1.41d3b43b744f48566cd3bb51015fp-96 },
        .{ .hi = 0x1.390f6844a0b83029d5524cap-3, .lo = 0x1.1c6d3f5f5fe70cb3a00ca5d7bfe4p-96 },
        .{ .hi = 0x1.3def0e6dfdf84ea10095aeap-3, .lo = 0x1.50447fcdd80256ad9e46d363cc54p-96 },
        .{ .hi = 0x1.42c7e7fe3fc01c5baa84ea8p-3, .lo = 0x1.e4a8511bcd6d383877c45cbc373cp-98 },
        .{ .hi = 0x1.479a07d3b641142ca3a5b06p-3, .lo = -0x1.e584ff398634ab49e12b7ceedcep-96 },
        .{ .hi = 0x1.4c65807e9333821962bd398p-3, .lo = -0x1.2b6fc634c578d7bb73b5f89170ecp-96 },
        .{ .hi = 0x1.512a644296c3cb096f47256p-3, .lo = -0x1.8eec43b9a26313b25a668455ed84p-100 },
        .{ .hi = 0x1.55e8c518b10f859bf037504p-3, .lo = 0x1.a44e1877778e042ca1fff8c09b4ap-97 },
        .{ .hi = 0x1.5aa0b4b0988f98f4b7b3558p-3, .lo = -0x1.1d36883ff38b16e03d088056210ap-97 },
        .{ .hi = 0x1.5f52447255c924e6e695998p-3, .lo = -0x1.c9a1067c1f62163817e106dfbadep-101 },
        .{ .hi = 0x1.63fd857fc49baa7c0cd1066p-3, .lo = 0x1.a97de9494385e4a74667a088fbep-96 },
        .{ .hi = 0x1.68a288b60b7fc2b622430e6p-3, .lo = -0x1.7f3541588eff7daf9b0d0b8edd73p-96 },
        .{ .hi = 0x1.6d415eaf0906a9ea9d132f2p-3, .lo = -0x1.c702d2f5f2dcb25f0ecb994061c2p-97 },
        .{ .hi = 0x1.71da17c2b7e7fea4e079fe8p-3, .lo = -0x1.47eea7b2200159b5df6c05428549p-101 },
        .{ .hi = 0x1.766cc40889e84a226ff02fp-3, .lo = 0x1.917c1270bcfda77828cfd78f80afp-96 },
        .{ .hi = 0x1.7af97358b9e03ccb5c4489p-3, .lo = 0x1.ef43da07a075a4654b00b1f15f38p-96 },
        .{ .hi = 0x1.7f80354d9529f92f3616978p-3, .lo = -0x1.1c614b89ec92ffd84e0260d9ac71p-98 },
        .{ .hi = 0x1.84011944bcb752c5b993p-3, .lo = 0x1.f6379a5ccc5d7e42c37fc0437b5p-97 },
        .{ .hi = 0x1.887c2e605e1189c603f25d8p-3, .lo = 0x1.d8ba842be4586d7a18917d454d5fp-96 },
        .{ .hi = 0x1.8cf183886480c9b28b1f97ep-3, .lo = -0x1.13ada343a03c927e5baecee4f8dap-99 },
        .{ .hi = 0x1.9161276ba29783a4f607cb8p-3, .lo = -0x1.0402e057ffccff9044bfb591e807p-99 },
        .{ .hi = 0x1.95cb2880f45ba6eadb35486p-3, .lo = -0x1.f8c898594a4ac24fe8b17309e439p-96 },
        .{ .hi = 0x1.9a2f95085a45b927e60039p-3, .lo = 0x1.4200f26808586e68a1259f2f8f61p-98 },
        .{ .hi = 0x1.9e8e7b0c0d4be203de57e9ap-3, .lo = 0x1.a25d8740fd57f8c816266d1e17dap-98 },
        .{ .hi = 0x1.a2e7e8618c2d24882a4f9dep-3, .lo = 0x1.42b5fc88aa72db35ef217dfcc3fdp-97 },
        .{ .hi = 0x1.a73beaaaa22f38e04a37a22p-3, .lo = -0x1.a404bfe25ad67d315d8ab0dd67dbp-98 },
        .{ .hi = 0x1.ab8a8f56677fc365b0e5a08p-3, .lo = -0x1.15fd6e4c7bf48b677423f326e48ep-97 },
        .{ .hi = 0x1.afd3e3a23b6800f54642bb8p-3, .lo = -0x1.b0ab128ccd50f29b68fe6c9d1a42p-97 },
        .{ .hi = 0x1.b417f49ab8806bc9543817ap-3, .lo = 0x1.19354df84685acaf5f6acdc7ba41p-103 },
        .{ .hi = 0x1.b856cf1ca31056c3f6e26b8p-3, .lo = -0x1.9bad52e70273c0d62c3c1c56dffdp-96 },
        .{ .hi = 0x1.bc907fd5d1c4069339ff22ep-3, .lo = 0x1.55d384963ee98afc3078d3044017p-99 },
        .{ .hi = 0x1.c0c5134610e267bfa808f84p-3, .lo = 0x1.c81a730b20e77706d7095ec3d1ccp-97 },
        .{ .hi = 0x1.c4f495c0002a25ee9a870fep-3, .lo = -0x1.710df049128528f24a199807e407p-96 },
        .{ .hi = 0x1.c91f1369eb7c9ad8af7db3ap-3, .lo = 0x1.d153b6ad5083f5c067d53c69945dp-98 },
        .{ .hi = 0x1.cd44983e9e7bca1ed1e0c96p-3, .lo = 0x1.f1cd40bfbf47786541d7b909c8f8p-96 },
        .{ .hi = 0x1.d165300e333f69c028a3c44p-3, .lo = 0x1.7ca1f33a3faaee8e8cefe3aea277p-96 },
        .{ .hi = 0x1.d580e67edc43ccfa0daf23p-3, .lo = 0x1.391a32168db31160772fa321faeap-98 },
        .{ .hi = 0x1.d997c70da9b46857c60d08ep-3, .lo = 0x1.fafa52e38c0cfb91035ed26ec006p-98 },
        .{ .hi = 0x1.dda9dd0f4a329136847dd6ap-3, .lo = -0x1.772bf9a47989f5ddf97c379a24afp-96 },
        .{ .hi = 0x1.e1b733b0c7381094f8c216p-3, .lo = -0x1.1f64198a27f7644abf7fc0a11cfep-100 },
        .{ .hi = 0x1.e5bfd5f83d342043025796cp-3, .lo = 0x1.0f7719e267ad0293ec177087d649p-96 },
        .{ .hi = 0x1.e9c3cec58f8072098058f2cp-3, .lo = -0x1.537f665ddb305fd1967d9b8f6ac6p-96 },
        .{ .hi = 0x1.edc328d3184af1cba1b7464p-3, .lo = 0x1.b5dd2a5b9ae43857a37357d141f9p-96 },
        .{ .hi = 0x1.f1bdeeb654900d96cd34f7ep-3, .lo = -0x1.5fbba7898678670262ca8d3b731dp-102 },
        .{ .hi = 0x1.f5b42ae08c4070bc91804dep-3, .lo = -0x1.269e7fd5a55d2611a378ea9f1f31p-96 },
        .{ .hi = 0x1.f9a5e79f76ac491748bb64p-3, .lo = 0x1.46fc084fcb735851dd511e9ab6fap-99 },
        .{ .hi = 0x1.fd932f1ddb4d5f2e278b32cp-3, .lo = 0x1.c1003738b4ba05532a543b03e65p-97 },
        .{ .hi = 0x1.00be05b217844161e1a46ep-2, .lo = -0x1.be23b7ac1afb629cbb08936bc5dep-95 },
        .{ .hi = 0x1.02b0432c96ff0694c1c8d1p-2, .lo = 0x1.f6d14c8185bcf5442230175deda1p-96 },
        .{ .hi = 0x1.04a054e13900409a780a5b4p-2, .lo = -0x1.9810c5ed46a87b19c7c9d5bf41bfp-96 },
        .{ .hi = 0x1.068e3fa282e3ced3324274cp-2, .lo = 0x1.3a690ff4678a2d9a79a0e0bdd6a7p-95 },
        .{ .hi = 0x1.087a0832fa7ac4f9ab7853ep-2, .lo = 0x1.2ddeb9fa1dc34acc69dec5d2cb6ap-95 },
        .{ .hi = 0x1.0a63b3456c818f3ddb757ccp-2, .lo = 0x1.6ef2b857afbb6f30837738a29904p-95 },
        .{ .hi = 0x1.0c4b457d3193d3ffa651b8cp-2, .lo = -0x1.dc7e54b397fe0d0f8c8247110d37p-96 },
        .{ .hi = 0x1.0e30c36e71a7f53a9ae38e2p-2, .lo = -0x1.3f13cdf87e3cc71194baf08e01e9p-96 },
        .{ .hi = 0x1.1014319e661bc87f6e8c7fep-2, .lo = -0x1.fe99c60b51d6330f43bbbc87a621p-98 },
        .{ .hi = 0x1.11f594839a5bd3aec4ea7c4p-2, .lo = 0x1.9056df9f7cd47dc0aaa4fe49396p-96 },
        .{ .hi = 0x1.13d4f0862b2e167244a4e9ap-2, .lo = -0x1.0ed925baabe232dd347a41af4531p-95 },
        .{ .hi = 0x1.15b24a0004a924955aced3ap-2, .lo = 0x1.b35a78737b603a1d2f1f8bbe02e1p-97 },
        .{ .hi = 0x1.178da53d1ee013c7b3e96d2p-2, .lo = 0x1.be3f91ccdbf2e4a1fce01e0b4bfcp-98 },
        .{ .hi = 0x1.1967067bb94b7feaa558f3p-2, .lo = -0x1.c81f29428e7ae8a88756077d6325p-95 },
        .{ .hi = 0x1.1b3e71ec94f7abbbb3324ap-2, .lo = 0x1.ff44989b978452e339c2f0e43dbep-95 },
        .{ .hi = 0x1.1d13ebb32d7f886517823d2p-2, .lo = 0x1.e1e49f3439556b32d3b34bb08983p-98 },
        .{ .hi = 0x1.1ee777e5f0dc35268e3c038p-2, .lo = 0x1.aa7013c4f4b0d97fb679b12eddecp-97 },
        .{ .hi = 0x1.20b91a8e761050d250ea2a8p-2, .lo = -0x1.0fb80164cc712614d5d1d7e782aep-99 },
        .{ .hi = 0x1.2288d7a9b2b641328381702p-2, .lo = 0x1.33609720003db913423ea05dd419p-96 },
        .{ .hi = 0x1.2456b3282f78608173a4448p-2, .lo = 0x1.2a9848b9f2603f28707cae68c496p-96 },
        .{ .hi = 0x1.2622b0ee3b79cee2bf4fc8ep-2, .lo = 0x1.2cc1e1efee6e9f2b64a0065118dbp-95 },
        .{ .hi = 0x1.27ecd4d41eb6752d3061152p-2, .lo = -0x1.27facfd964e8adddb3b8eaa2b26fp-95 },
        .{ .hi = 0x1.29b522a64b609745e857b1ep-2, .lo = 0x1.cc4f716351f4a4884ecf39351157p-96 },
        .{ .hi = 0x1.2b7b9e258e4226bf0485faap-2, .lo = 0x1.f4f31c49c1c1d7cfa213f77ce21p-95 },
        .{ .hi = 0x1.2d404b073e27da5069cad6ep-2, .lo = 0x1.7b2c22fa544845510f38e284490fp-95 },
        .{ .hi = 0x1.2f032cf56a5be40baedb9a4p-2, .lo = 0x1.3c375671456fd078abda0409a5c6p-95 },
        .{ .hi = 0x1.30c4478f0835f6cf717bf68p-2, .lo = -0x1.a5110af9dfae1fed23bb42c31e6dp-95 },
        .{ .hi = 0x1.32839e681fc6236e91f3dacp-2, .lo = 0x1.2bae43ce02a91a850fe7572c9b34p-95 },
        .{ .hi = 0x1.34413509f79fef311f12b36p-2, .lo = -0x1.fa41b743eca9679d5e17065b3ad1p-96 },
    },
};
