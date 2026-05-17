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
      label: 'Install',
      items: [
        {
          type: 'doc',
          id: 'installation/docker',
          label: 'Docker',
        },
        {
          type: 'doc',
          id: 'installation/google-cloud-docker',
          label: 'Google Cloud',
        },
        {
          type: 'doc',
          id: 'installation/apptainer',
          label: 'Apptainer / Singularity',
        },
        {
          type: 'doc',
          id: 'installation/non-containerized',
          label: 'From Source',
        },
      ],
    },
    {
      type: 'category',
      label: 'Run',
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
          label: 'WES Example',
        },
        {
          type: 'doc',
          id: 'usage/end-to-end-example-mit',
          label: 'mtDNA Example',
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
      label: 'Reproducibility',
      items: [
        {
          type: 'doc',
          id: 'usage/run-comparison',
          label: 'Run Comparison',
        },
        {
          type: 'doc',
          id: 'usage/integration-tests',
          label: 'Integration Tests',
        },
        {
          type: 'doc',
          id: 'usage/resource-validation',
          label: 'Resource Validation',
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
      label: 'Extend',
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
          type: 'doc',
          id: 'technical-details/adding-resources',
          label: 'Adding Resources',
        },
        {
          type: 'doc',
          id: 'technical-details/resource-bundle-v1',
          label: 'Resource Bundle v1',
        },
      ],
    },
    {
      type: 'category',
      label: 'Reference',
      items: [
        {
          type: 'doc',
          id: 'help/configuration-reference',
          label: 'Configuration',
        },
        {
          type: 'doc',
          id: 'help/outputs',
          label: 'Outputs',
        },
        {
          type: 'doc',
          id: 'help/naming-conventions',
          label: 'Naming Conventions',
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
