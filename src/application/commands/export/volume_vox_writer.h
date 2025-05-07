#ifndef CUBIQUITY_VOLUME_VOX_WRITER_H
#define CUBIQUITY_VOLUME_VOX_WRITER_H

#include "vox_writer/vox_writer.h"

#include "base/metadata.h"
#include "base/types.h"

#include "storage.h"

class volume_vox_writer : public vox_writer
{
public:
	volume_vox_writer(Cubiquity::Volume& vol, const Metadata& metadata);

protected:
	box  bounds() override;
	u8 voxel(const vec3i& position) override;

	void    on_progress(int done, int total) override;

private:
	Cubiquity::Volume& m_vol;
	const Metadata& m_metadata;
};

#endif // CUBIQUITY_VOLUME_VOX_WRITER_H
