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

fn __exp10h(x: f16) callconv(.c) f16 {
    _ = x;
    return 0;
}

fn exp10f(x: f32) callconv(.c) f32 {
    _ = x;
    return 0;
}

fn exp10(x: f64) callconv(.c) f64 {
    _ = x;
    return 0;
}

fn __exp10x(x: f80) callconv(.c) f80 {
    _ = x;
    return 0;
}

fn exp10q(x: f128) callconv(.c) f128 {
    _ = x;
    return 0;
}

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

