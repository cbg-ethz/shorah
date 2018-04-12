#!/bin/sh

install "${MESON_BUILD_ROOT}/cli.py" \
    "${MESON_BUILD_ROOT}/__main__.py" \
    "${MESON_BUILD_ROOT}/shotgun.py" \
    "${MESON_BUILD_ROOT}/amplicon.py" \
    "${MESON_BUILD_ROOT}/shorah_snv.py" \
    "${MESON_SOURCE_ROOT}/src/shorah"

install "${MESON_BUILD_ROOT}/src/cpp/b2w" \
    "${MESON_BUILD_ROOT}/src/cpp/diri_sampler" \
    "${MESON_BUILD_ROOT}/src/cpp/fil" \
    "${MESON_SOURCE_ROOT}/src/shorah/bin"
