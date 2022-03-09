#!/usr/bin/env bash

# don't want dependencies
cargo doc --no-deps
# remove old docs
rm -rf ./docs
# magic..?
echo "<meta http-equiv=\"refresh\" content=\"0; url=build_wheel\">" > target/doc/index.html
# copy to docs
cp -r target/doc ./docs