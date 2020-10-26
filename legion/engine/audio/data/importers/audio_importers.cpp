#include <audio/data/importers/audio_importers.hpp>
#if !defined(DOXY_EXCLUDE)
#define MINIMP3_IMPLEMENTATION
#include <minimp3.h>
#include <minimp3_ex.h>
#endif
#include <audio/systems/audiosystem.hpp>

namespace legion::audio
{
    common::result_decay_more<audio_segment, fs_error> mp3_audio_loader::load(const fs::basic_resource& resource, audio_import_settings&& settings)
    {
        using common::Err, common::Ok;
        using decay = common::result_decay_more<audio_segment, fs_error>;

        mp3dec_map_info_t map_info;
        map_info.buffer = resource.data();
        map_info.size = resource.size();

        mp3dec_t mp3dec;
        mp3dec_file_info_t fileInfo;

        if (mp3dec_load_mapinfo(&mp3dec, &map_info, &fileInfo, NULL, NULL))
        {
            return decay(Err(legion_fs_error("Failed to load audio file")));
        }

        audio_segment as(
            new byte[fileInfo.samples * 2], // fileInfo.samples is int16, therefore byte requires twice as much
            0,
            fileInfo.samples,
            fileInfo.channels,
            fileInfo.hz,
            fileInfo.layer,
            fileInfo.avg_bitrate_kbps
        );
        memmove(as.getData(), fileInfo.buffer, as.samples * sizeof(int16));
        free(fileInfo.buffer);

        async::readwrite_guard guard(AudioSystem::contextLock);
        alcMakeContextCurrent(AudioSystem::alcContext);
        //Generate openal buffer
        alGenBuffers((ALuint)1, &as.audioBufferId);
        ALenum format = AL_FORMAT_MONO16;
        if (as.channels == 2) format = AL_FORMAT_STEREO16;
        alBufferData(as.audioBufferId, format, as.getData(), as.samples * sizeof(int16), as.sampleRate);

        alcMakeContextCurrent(nullptr);

        return decay(Ok(as));
    }

    common::result_decay_more<audio_segment, fs_error> wav_audio_loader::load(const fs::basic_resource& resource, audio_import_settings&& settings)
    {
        using common::Err, common::Ok;
        using decay = common::result_decay_more<audio_segment, fs_error>;

        RIFF_Header header;
        WAVE_Data waveData;

        memcpy(&header, resource.data() , sizeof(header)); // Copy header data into the header struct

        // Check if the loaded file has the correct header

        if (header.chunckId[0] != 'R' ||
            header.chunckId[1] != 'I' ||
            header.chunckId[2] != 'F' ||
            header.chunckId[3] != 'F')
        {
            return decay(Err(legion_fs_error("WAV File invalid header, exptected RIFF")));
        }

        if (header.format[0] != 'W' ||
            header.format[1] != 'A' ||
            header.format[2] != 'V' ||
            header.format[3] != 'E')
        {
            return decay(Err(legion_fs_error("Loaded File is not of type WAV")));
        }

        if (header.wave_format.subChunckId[0] != 'f' ||
            header.wave_format.subChunckId[1] != 'm' ||
            header.wave_format.subChunckId[2] != 't' ||
            header.wave_format.subChunckId[3] != ' ')
        {
            return decay(Err(legion_fs_error("WAV File sub chunck id was not (fmt )")));
        }

        memcpy(&waveData, resource.data() + sizeof(header), sizeof(waveData));
        if (waveData.subChunckId[0] != 'd' ||
            waveData.subChunckId[1] != 'a' ||
            waveData.subChunckId[2] != 't' ||
            waveData.subChunckId[3] != 'a')
        {
            return decay(Err(legion_fs_error("WAV File sample data does not start with word (data)")));
        }

        const unsigned long sampleDataSize = resource.size() - sizeof(header) - sizeof(waveData);
        byte* audioData = new byte[sampleDataSize];
        memcpy(audioData, resource.data() + sizeof(header) + sizeof(waveData), sampleDataSize);
        
        audio_segment as(
            new byte[sampleDataSize],
            0,
            -1, // Sample count, unknown for wav
            header.wave_format.channels,
            (int)header.wave_format.sampleRate,
            -1, // Layer, does not exist in wav
            -1 // avg_biterate_kbps, unknown for wav
        );
        memmove(as.getData(), audioData, sampleDataSize);
        delete[] audioData; // Data has been moved, thus old data can be deleted

        async::readwrite_guard guard(AudioSystem::contextLock);
        alcMakeContextCurrent(AudioSystem::alcContext);
        //Generate openal buffer
        alGenBuffers((ALuint)1, &as.audioBufferId);
        ALenum format = AL_FORMAT_MONO16;
        if (as.channels == 1)
        {
            if (header.wave_format.bitsPerSample == 8) format = AL_FORMAT_MONO8;
            if (header.wave_format.bitsPerSample == 16) format = AL_FORMAT_MONO16;
            if (header.wave_format.bitsPerSample == 32) format = AL_FORMAT_MONO_FLOAT32;
        }
        else if (as.channels == 2)
        {
            if (header.wave_format.bitsPerSample == 8) format = AL_FORMAT_STEREO8;
            if (header.wave_format.bitsPerSample == 16) format = AL_FORMAT_STEREO16;
            if (header.wave_format.bitsPerSample == 32) format = AL_FORMAT_STEREO_FLOAT32;
        }
        alBufferData(as.audioBufferId, format, as.getData(), sampleDataSize, as.sampleRate);

        alcMakeContextCurrent(nullptr);

        return decay(Ok(as));
    }
}
