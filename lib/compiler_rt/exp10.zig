// f32 and f64 functions ported from
// musl/src/math/exp10f.c and
// musl/src/math/exp10.c

const std = @import("std");
const math = std.math;

const common = @import("common.zig");

comptime {
    @export(&__exp10h, .{ .name = "__exp10h", .linkage = common.linkage, .visibility = common.visibility });
    @export(&exp10f, .{ .name = "exp10f", .linkage = common.linkage, .visibility = common.visibility });
    @export(&exp10, .{ .name = "exp10", .linkage = common.linkage, .visibility = common.visibility });
    @export(&__exp10x, .{ .name = "__exp10x", .linkage = common.linkage, .visibility = common.visibility });
    if (common.want_ppc_abi) {
        @export(&exp10q, .{ .name = "exp10f128", .linkage = common.linkage, .visibility = common.visibility });
    }
    @export(&exp10q, .{ .name = "exp10q", .linkage = common.linkage, .visibility = common.visibility });
    @export(&exp10l, .{ .name = "exp10l", .linkage = common.linkage, .visibility = common.visibility });
}

const log2_10 = 3.3219280948873623478703194294893901758648313930246;

fn __exp10h(x: f16) callconv(.c) f16 {
    // TODO: more efficient implementation
    return @floatCast(exp10f(x));
}

fn exp10f(x: f32) callconv(.c) f32 {
    const p10 = [_]f32{ //
        1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1,
        1e1,  1e2,  1e3,  1e4,  1e5,  1e6,  1e7,
    };

    const ny = math.modf(x);
    const n = ny.ipart;
    const y = ny.fpart;
    const m: i32 = @intFromFloat(n);
    const i: u32 = @bitCast(n);
    // fabsf(n) < 8 without raising invalid on nan
    if (((i >> 23) & 0xff) < 0x7f + 3) {
        const p = p10[@intCast(m + 7)];
        if (y == 0) return p;
        return @exp2(log2_10 * y) * p;
    }

    const z: f64 = x;
    return @floatCast(@exp2(log2_10 * z));
}

fn exp10(x: f64) callconv(.c) f64 {
    const p10 = [_]f32{ //
        1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7,
        1e-6,  1e-5,  1e-4,  1e-3,  1e-2,  1e-1,  1,    1e1,  1e2,
        1e3,   1e4,   1e5,   1e6,   1e7,   1e8,   1e9,  1e10, 1e11,
        1e12,  1e13,  1e14,  1e15,
    };

    const ny = math.modf(x);
    const n = ny.ipart;
    const y = ny.fpart;
    const m: i64 = @intFromFloat(n);
    const i: u64 = @bitCast(n);

    // fabs(n) < 16 without raising invalid on nan
    if (((i >> 52) & 0x7ff) < 0x3ff + 4) {
        const p = p10[@intCast(m + 15)];
        if (y == 0) return p;
        return @exp2(log2_10 * y) * p;
    }
    return math.pow(f64, 10.0, x);
}

fn __exp10x(x: f80) callconv(.c) f80 {
    // TODO: more efficient implementation
    return @floatCast(exp10q(x));
}

pub const exp10q = @import("exp_f128.zig").exp10q;

fn exp10l(x: c_longdouble) callconv(.c) c_longdouble {
    switch (@typeInfo(c_longdouble).float.bits) {
        16 => return __exp10h(x),
        32 => return exp10f(x),
        64 => return exp10(x),
        80 => return __exp10x(x),
        128 => return exp10q(x),
        else => @compileError("unreachable"),
    }
}

const expect = std.testing.expect;
const expectRel = std.testing.expectApproxEqRel;

test exp10f {
    const eps = 50 * math.floatEps(f32);

    try expect(exp10f(0) == 1);

    try expectRel(1.584893193, exp10f(0.2), eps);
    try expectRel(7.81627805, exp10f(0.893), eps);
    try expectRel(31.6227766, exp10f(1.5), eps);
    try expectRel(2.8183828546771552e+37, exp10f(37.45), eps);
    try expectRel(0.1, exp10f(-1), eps);
}

test exp10 {
    const eps = 50 * math.floatEps(f32);

    try expect(exp10(0) == 1);

    try expectRel(1.584893193, exp10(0.2), eps);
    try expectRel(7.81627805, exp10(0.893), eps);
    try expectRel(31.6227766, exp10(1.5), eps);
    try expectRel(2.8183828546771552e+37, exp10(37.45), eps);
    try expectRel(0.1, exp10(-1), eps);
}

test exp10q {
    const eps = 50 * math.floatEps(f32);

    try expect(exp10q(0) == 1);

    try expectRel(1.584893193, exp10q(0.2), eps);
    try expectRel(7.81627805, exp10q(0.893), eps);
    try expectRel(31.6227766, exp10q(1.5), eps);
    try expectRel(2.8183828546771552e+37, exp10q(37.45), eps);
    try expectRel(0.1, exp10q(-1), eps);
}
