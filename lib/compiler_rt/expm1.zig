const common = @import("common.zig");

comptime {
    @export(&__expm1h, .{ .name = "__expm1h", .linkage = common.linkage, .visibility = common.visibility });
    @export(&expm1f, .{ .name = "expm1f", .linkage = common.linkage, .visibility = common.visibility });
    @export(&expm1, .{ .name = "expm1", .linkage = common.linkage, .visibility = common.visibility });
    @export(&__expm1x, .{ .name = "__expm1x", .linkage = common.linkage, .visibility = common.visibility });
    if (common.want_ppc_abi) {
        @export(&expm1q, .{ .name = "expm1f128", .linkage = common.linkage, .visibility = common.visibility });
    }
    @export(&expm1q, .{ .name = "expm1q", .linkage = common.linkage, .visibility = common.visibility });
    @export(&expm1l, .{ .name = "expm1l", .linkage = common.linkage, .visibility = common.visibility });

    // should not be needed
    @export(&__expm1h, .{ .name = "llvm.expm1.f16", .linkage = common.linkage, .visibility = common.visibility });
    @export(&expm1f, .{ .name = "llvm.expm1.f32", .linkage = common.linkage, .visibility = common.visibility });
    @export(&expm1, .{ .name = "llvm.expm1.f64", .linkage = common.linkage, .visibility = common.visibility });
    @export(&expm1q, .{ .name = "llvm.expm1.f128", .linkage = common.linkage, .visibility = common.visibility });
}

pub fn __expm1h(a: f16) callconv(.c) f16 {
    _ = a;
    return 0;
}
pub fn expm1f(a: f32) callconv(.c) f32 {
    _ = a;
    return 0;
}
pub fn expm1(a: f64) callconv(.c) f64 {
    _ = a;
    return 0;
}

pub fn __expm1x(a: f80) callconv(.c) f80 {
    return @floatCast(expm1q(a));
}

pub fn expm1q(a: f128) callconv(.c) f128 {
    _ = a;
    return 0;
}

pub fn expm1l(a: c_longdouble) callconv(.c) c_longdouble {
    switch (@typeInfo(c_longdouble).float.bits) {
        16 => return __expm1h(a),
        32 => return expm1f(a),
        64 => return expm1(a),
        80 => return __expm1x(a),
        128 => return expm1q(a),
        else => @compileError("unreachable"),
    }
}
