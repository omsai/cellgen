# Build the `fast_delaunay` executable on GNU/Linux.

package = fast-poisson-disk
url = https://github.com/thouis/$(package)/archive/d4834e72ea47837dcdc557c4ab4014498e2870f2.tar.gz
tarball = $(package)-$(lastword $(subst /, ,$(url)))
package_dir = $(tarball:.tar.gz=)
work_dir = $(package)

$(work_dir)/fast_delaunay : $(work_dir)/Makefile
	make -C $(@D)

$(work_dir)/Makefile : $(tarball) $(package)_make_linux.patch
	tar -xvpf $<
	touch -r $(package_dir) $<
	ln -s $(package_dir) $(@D)
	patch -d $(@D) -p1 < $(lastword $^)

$(tarball) :
	wget $(url) -O $@
