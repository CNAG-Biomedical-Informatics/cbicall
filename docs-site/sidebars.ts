import type {SidebarsConfig} from '@docusaurus/plugin-content-docs';

const sidebars: SidebarsConfig = {
  docsSidebar: [
    {
      type: 'doc',
      id: 'overview',
      label: 'Overview',
    },
    {
      type: 'category',
      label: 'Installation',
      items: [
        {
          type: 'doc',
          id: 'installation/docker',
          label: 'Docker',
        },
        {
          type: 'doc',
          id: 'installation/google-cloud-docker',
          label: 'Google Cloud with Docker',
        },
        {
          type: 'doc',
          id: 'installation/apptainer',
          label: 'HPC with Apptainer / Singularity',
        },
        {
          type: 'doc',
          id: 'installation/non-containerized',
          label: 'Non-containerized',
        },
      ],
    },
    {
      type: 'category',
      label: 'Getting Started',
      items: [
        {
          type: 'doc',
          id: 'usage/choose-your-path',
          label: 'Choose Your Path',
        },
        {
          type: 'doc',
          id: 'usage/quickstart',
          label: 'Quickstart',
        },
        {
          type: 'doc',
          id: 'usage/end-to-end-example-wes',
          label: 'End-to-end Example: WES',
        },
        {
          type: 'doc',
          id: 'usage/end-to-end-example-mit',
          label: 'End-to-end Example: mtDNA',
        },
        {
          type: 'doc',
          id: 'usage/usage',
          label: 'General Usage',
        },
      ],
    },
    {
      type: 'category',
      label: 'Pipelines',
      items: [
        {
          type: 'doc',
          id: 'pipelines/wes-wgs-single',
          label: 'WES/WGS Single-Sample',
        },
        {
          type: 'doc',
          id: 'pipelines/wes-wgs-cohort',
          label: 'WES/WGS Cohort',
        },
        {
          type: 'doc',
          id: 'pipelines/mtdna',
          label: 'mtDNA',
        },
      ],
    },
    {
      type: 'category',
      label: 'Technical Details',
      items: [
        {
          type: 'doc',
          id: 'technical-details/architecture',
          label: 'Architecture',
        },
        {
          type: 'doc',
          id: 'technical-details/adding-a-pipeline',
          label: 'Adding a Pipeline',
        },
        {
          type: 'category',
          label: 'Resource Bundles',
          items: [
            {
              type: 'doc',
              id: 'technical-details/resource-bundle-v1',
              label: 'Resource Bundle v1',
            },
          ],
        },
      ],
    },
    {
      type: 'category',
      label: 'Help & Reference',
      items: [
        {
          type: 'doc',
          id: 'help/naming-conventions',
          label: 'Naming Conventions',
        },
        {
          type: 'doc',
          id: 'help/configuration-reference',
          label: 'Configuration Reference',
        },
        {
          type: 'doc',
          id: 'help/outputs',
          label: 'Outputs',
        },
        {
          type: 'doc',
          id: 'help/performance',
          label: 'Performance',
        },
        {
          type: 'doc',
          id: 'help/troubleshooting',
          label: 'Troubleshooting',
        },
        {
          type: 'doc',
          id: 'help/faq',
          label: 'FAQ',
        },
      ],
    },
    {
      type: 'category',
      label: 'About',
      items: [
        {
          type: 'doc',
          id: 'about/about',
          label: 'About',
        },
        {
          type: 'doc',
          id: 'about/citation',
          label: 'Citation',
        },
        {
          type: 'doc',
          id: 'about/disclaimer',
          label: 'Disclaimer',
        },
      ],
    },
  ],
};

export default sidebars;
